
    echo "使用方法: $0 <工作目录>"
    exit 1
fi

dir="$1"
summary_file="${dir}/summary_results.txt"

if [ ! -d "$dir" ]; then
    echo "错误: 目录 $dir 不存在!"
    exit 1
fi

FSC2=/EClab/Project/fanyihang/software/miniconda3/envs/fsc2_env/bin/fsc27093

> "$summary_file"

for i in {1..200}
do
    echo "开始运行第 $i 次重复..."
    cd "$dir" || { echo "无法进入目录 $dir"; exit 1; }
    
    run_dir="run$i"
    mkdir -p "$run_dir" || { echo "无法创建目录 $run_dir"; exit 1; }
    cd "$run_dir" || { echo "无法进入目录 $run_dir"; exit 1; }
    
    ln -sf "$dir"/Know*obs .
    ln -sf "$dir"/Know.est .
    ln -sf "$dir"/Know.tpl .
    
    echo "执行fastsimcoal2..."
    $FSC2 -e Know.est -t Know.tpl -m -n 100000 -L 40 -s 0 -M  -c 12 -0
    
    if [ ! -d "Know" ]; then
        echo "警告: 未找到Know目录，跳过R脚本执行"
        cd "$dir" || exit 1
        continue
    fi
    
    echo "执行R脚本分析..."
    cd Know || { echo "无法进入Know目录"; cd "$dir"; exit 1; }
    
    ln -sf ../Know*obs .
    ln -sf ../Know.est .
    ln -sf ../Know.tpl .
    ln -sf ../calai.new.R .
    
    Rscript calai.new.R Know
    
    aic_value=$(sed -n '2p' Know.AIC)
    
    bestlhoods_content=$(sed -n '2p' Know.bestlhoods)
    
    echo -e "run$i\t$aic_value\t$bestlhoods_content" >> "$summary_file"
    
    cd "$dir" || exit 1
    
    echo "第 $i 次重复完成"
    echo "----------------------------------------"
done

sed -i '1i run\tAIC\t'"$(head -n1 "${dir}/run1/Know/Know.bestlhoods")" "$summary_file"

echo "所有重复运行完成"
echo "汇总结果已保存至: $summary_file"
