import threading
def dump_thread():
    print
    for thr in threading.enumerate():
        print(thr)
dump_thread()
import pymol
pymol.pymol_launch=4
pymol.pymol_argv = [ 'pymol', '-qc'] #  Quiet / no GUI
from pymol import cmd
pymol.finish_launching()
dump_thread()

cmd.fetch('1TUP', async=False)
cmd.disable('all')
cmd.enable('1TUP')
cmd.bg_color('white')
cmd.hide('all')
cmd.show('cartoon')
#cmd.hide('cartoon', 'chain E+F')
#cmd.show('ribbon', 'chain E+F')
cmd.select('zinc', 'name zn')
cmd.show('sphere', 'zinc')
cmd.set('ray_trace_mode', 3)
cmd.png('1TUP.png', width=1980, height=1080, quiet=0, ray=1, prior=False)
dump_thread()

cmd.set('ray_trace_mode', 1)
cmd.png('TUP.png', width=1980, height=1080, quiet=0, ray=1, prior=False)
cmd.quit()
