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