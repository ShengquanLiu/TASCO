import os

for filename in os.listdir('.'):
    # 只处理包含 sim-4 的文件
    if 'sim-4' in filename:
        new_name = filename.replace('sim-4', 'sim-3')
        
        if filename != new_name:
            # 检查是否存在同名的 sim-3 文件
            if os.path.exists(new_name):
                print(f"⚠️ 冲突：{new_name} 已存在！跳过 {filename}，防止覆盖。")
                # 如果你确认要覆盖，可以注释掉上面的 print，解除下面两行的注释
                # os.remove(new_name)
                # os.rename(filename, new_name)
            else:
                os.rename(filename, new_name)
                print(f"成功：{filename} -> {new_name}")