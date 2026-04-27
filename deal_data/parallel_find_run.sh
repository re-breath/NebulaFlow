#!/usr/bin/env bash
# =============================================================================
# parallel_find_run.sh
# 在当前目录下递归查找指定文件名，进入其所在目录并行执行指定命令
#
# 用法：
#   parallel_find_run.sh <文件名模式> <命令> [选项]
#
# 示例：
#   parallel_find_run.sh "dump.xyz" "analyze_thermo_out"
#   parallel_find_run.sh "*.log"    "grep ERROR"
#   parallel_find_run.sh "dump.xyz" "analyze_thermo_out" --jobs 8
#   parallel_find_run.sh "dump.xyz" "analyze_thermo_out" --dry-run
#   parallel_find_run.sh "dump.xyz" "analyze_thermo_out" --search-dir /data
# =============================================================================

set -euo pipefail

# ── 颜色输出 ──────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

info()    { echo -e "${CYAN}[INFO]${RESET}  $*"; }
success() { echo -e "${GREEN}[OK]${RESET}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
error()   { echo -e "${RED}[ERR]${RESET}   $*" >&2; }

# ── 使用帮助 ──────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
${BOLD}用法：${RESET}
  $(basename "$0") <文件名模式> <命令> [选项]

${BOLD}参数：${RESET}
  <文件名模式>      要查找的文件名，支持 glob（如 "dump.xyz"、"*.log"）
  <命令>            在每个匹配文件所在目录执行的命令（命令中可用 {FILE} 代指文件名）

${BOLD}选项：${RESET}
  -j, --jobs N      并行任务数（默认：自动检测 CPU 核数）
  -d, --search-dir  搜索根目录（默认：当前目录）
  -s, --source      额外 source 的脚本（默认已包含 ~/.rebreath/rebreath-function）
  -n, --dry-run     只打印将要执行的命令，不实际运行
  -q, --quiet       不打印每个任务的输出（只显示进度）
  -h, --help        显示此帮助

${BOLD}示例：${RESET}
  $(basename "$0") "dump.xyz"  "analyze_thermo_out"
  $(basename "$0") "dump.xyz"  "analyze_thermo_out" --jobs 8
  $(basename "$0") "*.lammps" "lmp -in {FILE}"       # {FILE} 会被替换为文件名
  $(basename "$0") "dump.xyz"  "analyze_thermo_out" --dry-run

# 手动指定并发数
$(basename "$0") "dump.xyz" "analyze_thermo_out" --jobs 8

# 先预览不执行（dry-run）
$(basename "$0") "dump.xyz" "analyze_thermo_out" --dry-run

# 从指定目录开始搜索
$(basename "$0") "dump.xyz" "analyze_thermo_out" --search-dir /data/sim

# {FILE} 占位符：命令中需要引用文件名时
$(basename "$0") "*.lammps" "lmp -in {FILE}" 

# 静默模式（不打印命令输出，只看进度）
$(basename "$0") "dump.xyz" "analyze_thermo_out" --quiet

EOF
    exit 0
}

# ── 默认参数 ──────────────────────────────────────────────────────────────────
PATTERN=""
COMMAND=""
SEARCH_DIR="$(pwd)"
JOBS=0            # 0 = 自动
DRY_RUN=false
QUIET=false
EXTRA_SOURCE=""

# ── 解析参数 ──────────────────────────────────────────────────────────────────
[[ $# -lt 2 ]] && { usage; }

PATTERN="$1"; shift
COMMAND="$1"; shift

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)        JOBS="$2";        shift 2 ;;
        -d|--search-dir)  SEARCH_DIR="$2";  shift 2 ;;
        -s|--source)      EXTRA_SOURCE="$2";shift 2 ;;
        -n|--dry-run)     DRY_RUN=true;     shift   ;;
        -q|--quiet)       QUIET=true;       shift   ;;
        -h|--help)        usage ;;
        *) error "未知选项：$1"; usage ;;
    esac
done

# ── 检测 CPU 核数 ─────────────────────────────────────────────────────────────
detect_cores() {
    local cores
    # 优先使用可调度核数，其次物理核数，最后兜底为 4
    if command -v nproc &>/dev/null; then
        cores=$(nproc --all 2>/dev/null || nproc)
    elif [[ -f /proc/cpuinfo ]]; then
        cores=$(grep -c ^processor /proc/cpuinfo)
    elif command -v sysctl &>/dev/null; then
        cores=$(sysctl -n hw.logicalcpu 2>/dev/null || echo 4)
    else
        cores=4
    fi
    echo "$cores"
}

CPU_CORES=$(detect_cores)
[[ "$JOBS" -le 0 ]] && JOBS="$CPU_CORES"

# ── 查找文件 ──────────────────────────────────────────────────────────────────
info "搜索目录：${BOLD}${SEARCH_DIR}${RESET}"
info "文件模式：${BOLD}${PATTERN}${RESET}"
info "执行命令：${BOLD}${COMMAND}${RESET}"
info "CPU 核数：${CPU_CORES}  →  并行任务数：${BOLD}${JOBS}${RESET}"
echo ""

# 用 find 收集所有匹配路径
mapfile -t FILES < <(find "$SEARCH_DIR" -name "$PATTERN" -type f | sort)

TOTAL="${#FILES[@]}"
if [[ "$TOTAL" -eq 0 ]]; then
    warn "未找到任何匹配 \"${PATTERN}\" 的文件，退出。"
    exit 0
fi
info "共找到 ${BOLD}${TOTAL}${RESET} 个文件"
echo ""

$DRY_RUN && warn "【DRY-RUN 模式】以下命令不会实际执行："

# ── 单任务执行函数（供 export 给子 shell）────────────────────────────────────
run_one() {
    local filepath="$1"
    local cmd="$2"
    local quiet="$3"
    local dry="$4"
    local extra_src="$5"
    local idx="$6"
    local total="$7"

    local dir
    dir="$(dirname "$filepath")"
    local filename
    filename="$(basename "$filepath")"

    # {FILE} 占位符替换为实际文件名
    local real_cmd="${cmd//\{FILE\}/$filename}"

    if [[ "$dry" == "true" ]]; then
        echo "  [${idx}/${total}] cd ${dir} && ${real_cmd}"
        return 0
    fi

    # 每个任务在独立子 shell 中运行，source 函数库后执行
    local log
    log=$(bash --login -c "
        # source 个人函数库（忽略不存在的情况）
        [[ -f \"\$HOME/.rebreath/rebreath-function\" ]] && \
            source \"\$HOME/.rebreath/rebreath-function\" 2>/dev/null || true
        # source 额外脚本
        [[ -n '${extra_src}' && -f '${extra_src}' ]] && \
            source '${extra_src}' 2>/dev/null || true
        cd '${dir}' || exit 1
        ${real_cmd}
    " 2>&1)

    local status=$?
    if [[ "$quiet" != "true" ]]; then
        # 加前缀方便区分不同任务的输出
        echo "$log" | sed "s|^|  [${idx}/${total}] ${filename}: |"
    fi

    if [[ $status -eq 0 ]]; then
        echo "__SUCCESS__ [${idx}/${total}] ${dir}"
    else
        echo "__FAILED__  [${idx}/${total}] ${dir}"
    fi
}

export -f run_one

# ── 并行执行引擎 ──────────────────────────────────────────────────────────────
SUCCESS_COUNT=0
FAIL_COUNT=0
DONE_COUNT=0

# 用关联数组跟踪后台 PID → 索引
declare -A PID_MAP   # pid -> idx
declare -A IDX_MAP   # pid -> filepath

# 进度条打印
print_progress() {
    local done=$1 total=$2 label=$3
    local bar_len=35
    local filled=$(( bar_len * done / total ))
    local bar
    bar=$(printf '█%.0s' $(seq 1 $filled 2>/dev/null) 2>/dev/null || true)
    local empty
    empty=$(printf '░%.0s' $(seq 1 $((bar_len - filled)) 2>/dev/null) 2>/dev/null || true)
    local pct=$(( done * 100 / total ))
    printf "\r  [%s%s] %3d%%  %d/%d  %s" \
        "$bar" "$empty" "$pct" "$done" "$total" "$label"
}

# 临时目录存放各任务输出
TMPDIR_TASKS=$(mktemp -d)
trap 'rm -rf "$TMPDIR_TASKS"' EXIT

idx=0
for filepath in "${FILES[@]}"; do
    idx=$(( idx + 1 ))
    outfile="${TMPDIR_TASKS}/${idx}.out"

    if $DRY_RUN; then
        run_one "$filepath" "$COMMAND" "$QUIET" "true" "$EXTRA_SOURCE" "$idx" "$TOTAL"
        continue
    fi

    # 启动后台任务，输出重定向到临时文件
    run_one "$filepath" "$COMMAND" "$QUIET" "false" "$EXTRA_SOURCE" "$idx" "$TOTAL" \
        > "$outfile" 2>&1 &
    pid=$!
    PID_MAP[$pid]=$idx
    IDX_MAP[$pid]=$filepath

    # 当后台任务数达到 JOBS 时，等待最早完成的一个
    while [[ ${#PID_MAP[@]} -ge $JOBS ]]; do
        for pid in "${!PID_MAP[@]}"; do
            if ! kill -0 "$pid" 2>/dev/null; then
                # 任务结束，读取结果
                wait "$pid" 2>/dev/null || true
                out_content=$(cat "${TMPDIR_TASKS}/${PID_MAP[$pid]}.out" 2>/dev/null || true)

                if echo "$out_content" | grep -q "^__SUCCESS__"; then
                    (( SUCCESS_COUNT++ ))
                else
                    (( FAIL_COUNT++ ))
                fi
                $QUIET || echo "$out_content" | grep -v "^__" 

                (( DONE_COUNT++ ))
                print_progress "$DONE_COUNT" "$TOTAL" \
                    "$(basename "${IDX_MAP[$pid]}")"

                unset PID_MAP[$pid]
                unset IDX_MAP[$pid]
                break
            fi
        done
        sleep 0.1
    done
done

# 等待剩余任务完成
if ! $DRY_RUN; then
    for pid in "${!PID_MAP[@]}"; do
        wait "$pid" 2>/dev/null || true
        out_content=$(cat "${TMPDIR_TASKS}/${PID_MAP[$pid]}.out" 2>/dev/null || true)

        if echo "$out_content" | grep -q "^__SUCCESS__"; then
            (( SUCCESS_COUNT++ ))
        else
            (( FAIL_COUNT++ ))
        fi
        $QUIET || echo "$out_content" | grep -v "^__"

        (( DONE_COUNT++ ))
        print_progress "$DONE_COUNT" "$TOTAL" \
            "$(basename "${IDX_MAP[$pid]}")"
    done
fi

# ── 汇总 ──────────────────────────────────────────────────────────────────────
if ! $DRY_RUN; then
    echo -e "\n"
    echo -e "  ${BOLD}────────── 执行汇总 ──────────${RESET}"
    echo -e "  总计：${TOTAL}  ${GREEN}成功：${SUCCESS_COUNT}${RESET}  ${RED}失败：${FAIL_COUNT}${RESET}"
    [[ $FAIL_COUNT -gt 0 ]] && \
        warn "有失败任务，可加 --quiet 关闭静默或检查对应目录"
fi


