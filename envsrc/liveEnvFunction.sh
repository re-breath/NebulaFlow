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