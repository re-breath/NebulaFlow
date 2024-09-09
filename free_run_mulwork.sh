#!/bin/bash
#本脚本用来在空闲的gpu上运行任务，检查是否有gpu空闲，如果有则在该gpu上运行任务
#例子：work=’cd workdir;bash run.sh/gpumd‘
#注意：对于nep，gpumd等任务，想要使其不与外界进行交互静默运行可以使用nohup nep/gpumd 2>&1 &
work1='bash run_gpumd_5.sh'
#work2='cd 64atoms0_5fs;bash run_gpumd_5.sh'


task_commands=()
for var in $(compgen -A variable | grep '^work'); do
  task_commands+=("${!var}")
done
check_and_run() {
for work in "${task_commands[@]}"; do
  while true; do
    for gpu_id in $(nvidia-smi --query-gpu=index --format=csv,noheader,nounits); do
      mem_used=$(nvidia-smi --id=$gpu_id --query-gpu=memory.used --format=csv,noheader,nounits)
      if [ $mem_used -lt 200 ]; then
          export CUDA_VISIBLE_DEVICES=$gpu_id
          echo "Running task on GPU $gpu_id"
          eval $work
          break 2
      fi
      done
      echo "No free GPU found, waiting..."
      sleep 60
  done
  date '+%Y-%m-%d %H:%M:%S' >> run_train-file.log
  echo -e "执行 ${work} \n" >> run_train-file.log
done
}

check_and_run