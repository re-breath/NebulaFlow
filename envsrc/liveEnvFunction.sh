# 该函数库用来存放生活中常用的一些函数，其中包含了视频处理等功能


compress_all_gif() {
    # 压缩当前目录下所有gif文件，压缩后文件名为原文件名_compressed.gif
    # 使用方式为：compress_all_gif 3  压缩大于3MB的文件，压缩到3MB
    local compress_size=${1:-3}  # 默认压缩阈值为3MB
    local DIR=$PWD
    for file in "$DIR"/*.gif; do
        # 检查文件是否存在
        if [ ! -f "$file" ]; then
            continue
        fi

        # 获取文件大小（单位：字节）
        size=$(stat -c%s "$file")

        # 计算阈值（MB转字节）
        local threshold=$((compress_size * 1024 * 1024))

        # 检查文件大小是否大于指定阈值
        if [ "$size" -gt "$threshold" ]; then
            echo "Compressing $file (Size: $(($size / 1024 / 1024))MB > $compress_size MB)"
            # 使用ffmpeg压缩文件到接近指定大小
            ffmpeg -i "$file" -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -pix_fmt yuv420p -r 15 -b:v 1M -movflags +faststart -fs ${threshold} "${file%.gif}_compressed.gif"
        else
            echo "$file is already less than or equal to $compress_size MB (Size: $(($size / 1024 / 1024))MB)"
        fi
    done
}

convert_video_to_gif() {
    # 将视频文件转换为GIF文件
    # 使用方式为：convert_video_to_gif input_video.mp4 output.gif
    local input_video="$1"
    local videoname=${input_video%.*}
    local output_gif="${2:-${videoname}.gif}"

    # 检查输入文件是否存在
    if [ ! -f "$input_video" ]; then
        echo "Error: Input video file '$input_video' does not exist."
        return 1
    fi

    # 使用ffmpeg将视频转换为GIF
    ffmpeg -i "$input_video" -vf "fps=10,scale=320:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 "$output_gif"

    # 检查转换是否成功
    if [ $? -eq 0 ]; then
        echo "Video successfully converted to GIF: $output_gif"
    else
        echo "Error: Failed to convert video to GIF."
        return 1
    fi
}

compare_folders() {
    # 检查参数数量
    if [ "$#" -ne 2 ]; then
        echo "Usage: compare_folders <folder1> <folder2>"
        return 1
    fi

    # 获取文件夹路径
    folder1="$1"
    folder2="$2"

    # 检查文件夹是否存在
    if [ ! -d "$folder1" ]; then
        echo "Error: Folder '$folder1' does not exist."
        return 1
    fi

    if [ ! -d "$folder2" ]; then
        echo "Error: Folder '$folder2' does not exist."
        return 1
    fi

    # 获取两个文件夹中的文件列表
    files1=$(ls "$folder1")
    files2=$(ls "$folder2")

    # 检查同名文件
    for file in $files1; do
        if [ -f "$folder1/$file" ] && [ -f "$folder2/$file" ]; then
            # 比较同名文件
            if ! diff -q "$folder1/$file" "$folder2/$file" > /dev/null; then
                echo "Files differ: $file"
            else
                echo "Files are identical: $file"
            fi
        fi
    done

    # 检查 folder2 中独有的文件
    for file in $files2; do
        if [ -f "$folder2/$file" ] && [ ! -f "$folder1/$file" ]; then
            echo "File only in $folder2: $file"
        fi
    done

    # 检查 folder1 中独有的文件
    for file in $files1; do
        if [ -f "$folder1/$file" ] && [ ! -f "$folder2/$file" ]; then
            echo "File only in $folder1: $file"
        fi
    done
}

# 该函数用来将markdown文件转换为rich格式的文本，可在linux终端直接查看
catmd() {
    local mdfile="${1:-README.md}"

    if [[ ! -f "$mdfile" ]]; then
        echo "File not found: $mdfile" >&2
        return 1
    fi

    python3 - "$mdfile" <<'PY'
import sys
from rich.console import Console
from rich.markdown import Markdown

console = Console()
mdfile = sys.argv[1]

with open(mdfile, "r", encoding="utf-8") as f:
    console.print(Markdown(f.read()))
PY
}

# 该函数用来将多个png文件合并为一个pdf文件
# 依赖安装：pip install pillow
png2pdf() {
  # usage: png2pdf in.png [out.pdf]  OR  png2pdf *.png out.pdf
  python - <<'PY' "$@"
import sys, os
from PIL import Image

args = sys.argv[1:]
if not args:
    print("Usage: png2pdf <in.png|*.png> [out.pdf]")
    sys.exit(2)

# if last arg is pdf -> output
if len(args) >= 2 and args[-1].lower().endswith(".pdf"):
    out = args[-1]
    ins = args[:-1]
else:
    out = (os.path.splitext(args[0])[0] + ".pdf") if len(args)==1 else "merged.pdf"
    ins = args

imgs = []
for p in ins:
    im = Image.open(p)
    # handle transparency -> white background
    if im.mode in ("RGBA", "LA") or (im.mode == "P" and "transparency" in im.info):
        bg = Image.new("RGB", im.size, (255,255,255))
        bg.paste(im.convert("RGBA"), mask=im.convert("RGBA").split()[-1])
        im = bg
    else:
        im = im.convert("RGB")
    imgs.append(im)

first, rest = imgs[0], imgs[1:]
first.save(out, "PDF", save_all=True, append_images=rest)
print("Saved:", out)
PY
}


pdf2eps() {
  local input_dir="${1:-.}"
  local output_dir="${2:-./eps_output}"

  if ! command -v pdftops >/dev/null 2>&1; then
    echo "错误: 没有找到 pdftops，请先安装 poppler。"
    return 1
  fi

  mkdir -p "$output_dir" || return 1

  find "$input_dir" -type f -iname "*.pdf" | while IFS= read -r pdf; do
    local rel name out
    name="$(basename "$pdf" .pdf)"
    out="$output_dir/${name}.eps"

    echo "转换: $pdf -> $out"
    pdftops -eps "$pdf" "$out"
  done
}


pdfs_to_eps_pages() {
  local input_dir="${1:-.}"
  local output_dir="${2:-./eps_output}"

  if ! command -v pdfinfo >/dev/null 2>&1 || ! command -v pdftops >/dev/null 2>&1; then
    echo "错误: 需要 pdfinfo 和 pdftops，请先安装 poppler。"
    return 1
  fi

  mkdir -p "$output_dir" || return 1

  find "$input_dir" -type f -iname "*.pdf" | while IFS= read -r pdf; do
    local name pages i out
    name="$(basename "$pdf" .pdf)"
    pages="$(pdfinfo "$pdf" 2>/dev/null | awk '/^Pages:/ {print $2}')"

    if [ -z "$pages" ]; then
      echo "跳过: 无法读取页数 $pdf"
      continue
    fi

    for ((i=1; i<=pages; i++)); do
      out="$output_dir/${name}_p${i}.eps"
      echo "转换: $pdf 第 $i 页 -> $out"
      pdftops -eps -f "$i" -l "$i" "$pdf" "$out"
    done
  done
}

txt2mp3() {
  # usage:
  #   txt2mp3 in.txt
  #   txt2mp3 a.txt b.txt c.txt
  #   txt2mp3 *.txt
  #   txt2mp3 -j 4 *.txt
  #   txt2mp3 -v zh-CN-XiaoxiaoNeural *.txt
  #
  # output:
  #   in.txt -> in.mp3

  python - <<'PY' "$@"
import sys, os, asyncio, argparse
from pathlib import Path

try:
    import edge_tts
except ImportError:
    print("Missing dependency: pip install edge-tts", file=sys.stderr)
    sys.exit(1)

def parse_args():
    p = argparse.ArgumentParser(
        prog="txt2mp3",
        description="Convert txt files to mp3 in parallel using edge-tts."
    )
    p.add_argument("files", nargs="+", help="Input .txt files")
    p.add_argument("-j", "--jobs", type=int, default=4, help="Parallel jobs, default: 4")
    p.add_argument("-v", "--voice", default="zh-CN-XiaoxiaoNeural", help="TTS voice")
    p.add_argument("-r", "--rate", default="+0%", help="Speech rate, e.g. +10%, -10%")
    p.add_argument("-V", "--volume", default="+0%", help="Volume, e.g. +10%, -10%")
    return p.parse_args()

async def convert_one(path, voice, rate, volume, sem):
    async with sem:
        p = Path(path)

        if not p.exists():
            print(f"Skip missing file: {p}", file=sys.stderr)
            return False

        if p.is_dir():
            print(f"Skip directory: {p}", file=sys.stderr)
            return False

        out = p.with_suffix(".mp3")

        try:
            text = p.read_text(encoding="utf-8").strip()
        except UnicodeDecodeError:
            text = p.read_text(encoding="utf-8-sig").strip()

        if not text:
            print(f"Skip empty file: {p}", file=sys.stderr)
            return False

        try:
            communicate = edge_tts.Communicate(
                text=text,
                voice=voice,
                rate=rate,
                volume=volume,
            )
            await communicate.save(str(out))
            print(f"Saved: {out}")
            return True
        except Exception as e:
            print(f"Failed: {p} -> {e}", file=sys.stderr)
            return False

async def main():
    args = parse_args()

    if args.jobs < 1:
        print("--jobs must be >= 1", file=sys.stderr)
        sys.exit(2)

    sem = asyncio.Semaphore(args.jobs)

    tasks = [
        convert_one(f, args.voice, args.rate, args.volume, sem)
        for f in args.files
    ]

    results = await asyncio.gather(*tasks)
    ok = sum(1 for x in results if x)

    if ok == 0:
        sys.exit(1)

if __name__ == "__main__":
    asyncio.run(main())
PY
}

# 膬换为 MP3 格式
to_mp3() {
  if [ $# -lt 1 ]; then
    echo "用法: to_mp3 输入文件 [输出文件.mp3]"
    return 1
  fi

  local input="$1"
  local output="$2"

  if [ ! -f "$input" ]; then
    echo "文件不存在: $input"
    return 1
  fi

  if ! command -v ffmpeg >/dev/null 2>&1; then
    echo "未找到 ffmpeg，请先安装："
    echo "  macOS: brew install ffmpeg"
    echo "  Ubuntu/Debian: sudo apt install ffmpeg"
    return 1
  fi

  if [ -z "$output" ]; then
    local base="${input%.*}"
    output="${base}.converted.mp3"
  fi

  ffmpeg -y \
    -i "$input" \
    -vn \
    -map 0:a:0 \
    -c:a libmp3lame \
    -b:a 128k \
    -ar 44100 \
    -ac 2 \
    "$output"

  echo "已转换为真正的 MP3: $output"
}


to_flac_light() {
  if [ $# -lt 1 ]; then
    echo "用法: to_flac_light 输入文件 [输出文件.flac] [模式]"
    echo "模式:"
    echo "  voice   默认，轻量人声：单声道 22050 Hz"
    echo "  voice32 稍高音质人声：单声道 32000 Hz"
    echo "  keep    尽量保持原始采样率和声道，但文件会很大"
    return 1
  fi

  local input="$1"
  local output="$2"
  local mode="${3:-voice}"

  if [ ! -f "$input" ]; then
    echo "文件不存在: $input"
    return 1
  fi

  if ! command -v ffmpeg >/dev/null 2>&1; then
    echo "未找到 ffmpeg，请先安装："
    echo "  macOS: brew install ffmpeg"
    echo "  Ubuntu/Debian: sudo apt install ffmpeg"
    return 1
  fi

  if [ -z "$output" ]; then
    local base="${input%.*}"
    output="${base}.light.flac"
  fi

  case "$mode" in
    voice)
      ffmpeg -y \
        -i "$input" \
        -vn \
        -map 0:a:0 \
        -ac 1 \
        -ar 22050 \
        -c:a flac \
        -compression_level 8 \
        "$output"
      ;;
    voice32)
      ffmpeg -y \
        -i "$input" \
        -vn \
        -map 0:a:0 \
        -ac 1 \
        -ar 32000 \
        -c:a flac \
        -compression_level 8 \
        "$output"
      ;;
    keep)
      ffmpeg -y \
        -i "$input" \
        -vn \
        -map 0:a:0 \
        -c:a flac \
        -compression_level 8 \
        "$output"
      ;;
    *)
      echo "未知模式: $mode"
      echo "可用模式: voice, voice32, keep"
      return 1
      ;;
  esac

  echo "已转换为 FLAC: $output"
}


novel_txt_to_mp3() {
  if [ $# -lt 1 ]; then
    echo "用法: novel_txt_to_mp3 小说.txt [输出目录] [并行数]"
    echo "示例: novel_txt_to_mp3 book.txt ./mp3 4"
    return 1
  fi

  local input="$1"
  local outdir="${2:-./novel_mp3}"
  local jobs="${3:-4}"

  local voice="${TTS_VOICE:-zh-CN-XiaoxiaoNeural}"
  local rate="${TTS_RATE:--5%}"
  local pitch="${TTS_PITCH:-+8Hz}"
  local volume="${TTS_VOLUME:-+0%}"

  if [ ! -f "$input" ]; then
    echo "文件不存在: $input"
    return 1
  fi

  if ! command -v python3 >/dev/null 2>&1; then
    echo "未找到 python3"
    return 1
  fi

  if ! command -v edge-tts >/dev/null 2>&1; then
    echo "未找到 edge-tts，请先安装："
    echo "  python3 -m pip install -U edge-tts"
    return 1
  fi

  mkdir -p "$outdir"

  local tmpdir
  tmpdir="$(mktemp -d)"
  local chapter_dir="$tmpdir/chapters"
  local manifest="$tmpdir/manifest.tsv"
  mkdir -p "$chapter_dir"

  python3 - "$input" "$chapter_dir" "$outdir" "$manifest" <<'PY'
import sys
import re
from pathlib import Path

input_file = Path(sys.argv[1])
chapter_dir = Path(sys.argv[2])
outdir = Path(sys.argv[3])
manifest = Path(sys.argv[4])

raw = input_file.read_bytes()

text = None
for enc in ("utf-8-sig", "utf-8", "gb18030", "gbk", "big5"):
    try:
        text = raw.decode(enc)
        break
    except UnicodeDecodeError:
        pass

if text is None:
    raise SystemExit("无法识别 txt 编码，请先转成 UTF-8 或 GB18030")

text = text.replace("\r\n", "\n").replace("\r", "\n")

chapter_re = re.compile(
    r"^\s*("
    r"第[零一二三四五六七八九十百千万两〇0-9０-９]+[章节回卷集部篇].{0,80}"
    r"|卷[零一二三四五六七八九十百千万两〇0-9０-９]+.{0,80}"
    r"|Chapter\s+[0-9IVXLCDM]+.{0,80}"
    r"|CHAPTER\s+[0-9IVXLCDM]+.{0,80}"
    r")\s*$"
)

def safe_name(s: str, max_len: int = 50) -> str:
    s = re.sub(r"[\\/:*?\"<>|]", "_", s)
    s = re.sub(r"\s+", " ", s).strip()
    s = s.strip("._ ")
    return s[:max_len] or "untitled"

lines = text.split("\n")
chapters = []
title = "正文"
buf = []

for line in lines:
    stripped = line.strip()
    if chapter_re.match(stripped):
        if buf and "\n".join(buf).strip():
            chapters.append((title, "\n".join(buf).strip()))
        title = stripped
        buf = [stripped]
    else:
        buf.append(line)

if buf and "\n".join(buf).strip():
    chapters.append((title, "\n".join(buf).strip()))

if len(chapters) <= 1:
    body = text.strip()
    chunk_size = 3500
    chapters = []
    for i in range(0, len(body), chunk_size):
        title = f"第{i // chunk_size + 1}段"
        chapters.append((title, body[i:i + chunk_size]))

with manifest.open("w", encoding="utf-8") as mf:
    for idx, (title, body) in enumerate(chapters, 1):
        name = f"{idx:02d}.{safe_name(title)}"
        txt_path = chapter_dir / f"{name}.txt"
        mp3_path = outdir / f"{name}.mp3"

        txt_path.write_text(body, encoding="utf-8")
        mf.write(f"{txt_path}\t{mp3_path}\n")

print(f"已切分章节数: {len(chapters)}")
PY

  _tts_one_chapter() {
    local txt="$1"
    local mp3="$2"

    if [ -s "$mp3" ]; then
      echo "已存在，跳过: $(basename "$mp3")"
      return 0
    fi

    echo "开始生成: $(basename "$mp3")"

    edge-tts \
      --voice "$TTS_VOICE" \
      --rate="$TTS_RATE" \
      --pitch="$TTS_PITCH" \
      --volume="$TTS_VOLUME" \
      --file "$txt" \
      --write-media "$mp3"

    local status=$?

    if [ $status -ne 0 ] || [ ! -s "$mp3" ]; then
      echo "失败: $(basename "$mp3")"
      rm -f "$mp3"
      return 1
    fi

    echo "完成: $(basename "$mp3")"
    return 0
  }

  export -f _tts_one_chapter
  export TTS_VOICE="$voice"
  export TTS_RATE="$rate"
  export TTS_PITCH="$pitch"
  export TTS_VOLUME="$volume"

  local fail_log="$tmpdir/fail.log"
  : > "$fail_log"

  while IFS=$'\t' read -r txt mp3; do
    printf '%s\0%s\0' "$txt" "$mp3"
  done < "$manifest" |
    xargs -0 -n 2 -P "$jobs" bash -c '
      _tts_one_chapter "$1" "$2" || echo "$2" >> "$FAIL_LOG"
    ' _

  if [ -s "$fail_log" ]; then
    echo
    echo "有部分章节生成失败："
    cat "$fail_log"
    echo
    echo "你可以降低并行数后重试，例如："
    echo "  novel_txt_to_mp3 \"$input\" \"$outdir\" 2"
    rm -rf "$tmpdir"
    return 1
  fi

  echo
  echo "全部完成，输出目录：$outdir"
  echo "生成文件数量：$(find "$outdir" -maxdepth 1 -type f -name '*.mp3' | wc -l)"

  rm -rf "$tmpdir"
}

get_map_distance() {
    # 该函数可以查询两地之间的距离
    local AMAP_KEY="ebb8ae1f79a839931faafd8dcd8fdcfe"
    
    if [ $# -lt 1 ]; then
        echo "使用方法:"
        echo "  1) 查两地距离: mapdist <出发地> <目的地>"
        echo "  2) 查当前到某地: mapdist <目的地>"
        return 1
    fi

    local PARAM_1="$1"
    local PARAM_2="$2"
    local DEFAULT_CITY=""

    # 如果只传了一个参数，说明 PARAM_1 是目的地，出发地需要自动定位
    if [ -z "$PARAM_2" ]; then
        echo "📡 未检测到出发地，正在自动定位你当前的位置..."
    else
        echo "正在查询 (${DEFAULT_CITY}): [$PARAM_1] 到 [$PARAM_2] ..."
    fi

    python3 - "$AMAP_KEY" "$PARAM_1" "$PARAM_2" "$DEFAULT_CITY" << 'EOF'
import sys
import urllib.request
import urllib.parse
import json

key = sys.argv[1]
p1 = sys.argv[2]
p2 = sys.argv[3]
city_limit = sys.argv[4]

def get_current_ip_location():
    """通过高德 IP 定位 API 获取当前位置"""
    url = f"https://restapi.amap.com/v3/ip?key={key}"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            if data.get('status') == '1' and'rectangle' in data and data['rectangle']:
                # IP 返回的是一个矩形范围(如 "116.32,39.96;116.36,39.99")，我们取中心点
                rect = data['rectangle'].split(';')[0].split(',')
                return f"{rect[0]},{rect[1]}", f"当前位置({data.get('city', '未知城市')})"
    except Exception:
        pass
    return None, None

def get_geo(address):
    """地理编码：文字地址转经纬度"""
    encoded_addr = urllib.parse.quote(address)
    encoded_city = urllib.parse.quote(city_limit)
    url = f"https://restapi.amap.com/v3/geocode/geo?key={key}&address={encoded_addr}&city={encoded_city}"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            if data.get('status') == '1' and data.get('geocodes'):
                return data['geocodes'][0]['location'], data['geocodes'][0]['formatted_address']
    except Exception:
        pass
    return None, None

# 判断是单参数（自动定位）还是双参数
if not p2:
    dest_name = p1
    ori_loc, ori_format = get_current_ip_location()
    if not ori_loc:
        print("❌ 自动定位失败，请手动输入出发地。")
        sys.exit(1)
else:
    origin_name = p1
    dest_name = p2
    ori_loc, ori_format = get_geo(origin_name)
    if not ori_loc:
        print(f"❌ 无法解析出发地 [{origin_name}]，请确认该城市是否存在此地名。")
        sys.exit(1)

# 解析目的地
dest_loc, dest_format = get_geo(dest_name)
if not dest_loc:
    print(f"❌ 无法解析目的地 [{dest_name}]。")
    sys.exit(1)

# 计算驾车距离 (type=1)
dist_url = f"https://restapi.amap.com/v3/distance?key={key}&origins={ori_loc}&destination={dest_loc}&type=1"
try:
    req = urllib.request.Request(dist_url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req) as response:
        res_data = json.loads(response.read().decode('utf-8'))
        if res_data.get('status') == '1' and res_data.get('results'):
            result = res_data['results'][0]
            meters = float(result['distance'])
            seconds = float(result['duration'])
            
            print(f"\n🚗 【高德地图查询结果 (智能定位版)】")
            print(f"📍 出发地: {ori_format}")
            print(f"📍 目的地: {dest_name} ({dest_format})")
            print("--------------------------------")
            print(f"🛣️  驾车距离: {meters/1000:.2f} 公里")
            print(f"⏱️  预计耗时: {int(seconds//60)} 分钟 (不含严重拥堵)")
        else:
            print(f"❌ 距离计算失败。原因: {res_data.get('info', '未知错误')}")
except Exception as e:
    print(f"❌ 网络请求异常: {e}")
EOF
}
alias mapdist=get_map_distance
