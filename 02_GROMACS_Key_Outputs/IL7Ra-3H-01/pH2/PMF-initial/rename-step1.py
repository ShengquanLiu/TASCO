import os

# 获取当前文件夹下的所有文件
for filename in os.listdir('.'):
    # 仅处理包含 awh96 的文件
    if 'awh96' in filename:
        # 替换字符
        new_name = filename.replace('awh96', 'awh').replace('sim-3', 'sim-2')
        # 避免重命名错误
        if filename != new_name:
            try:
                os.rename(filename, new_name)
                print(f"成功: {filename} -> {new_name}")
            except FileExistsError:
                print(f"警告: 目标文件 {new_name} 已存在，跳过 {filename}")
            except Exception as e:
                print(f"错误: 重命名 {filename} 失败，原因: {e}")