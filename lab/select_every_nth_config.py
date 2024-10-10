# 该脚本使用来处理轨迹类型的文件，按照n取一来保存到一个新的文件中
import ase.io

nth = 10
filename = 'dump.xyz'
def select_every_nth_config(filename, nth):
    traj = ase.io.read(filename, index = ':')
    new_traj = traj[::nth]
    ase.io.write('new_traj.xyz', new_traj)

select_every_nth_config(filename, nth)