import pymol
from pymol import cmd
#pymol.pymol_argv = [ 'pymol', '-qc'] #  Quiet / no GUI
pymol.finish_launching()

cmd.fetch('1TUP', async=False)

cmd.disable('all')
cmd.enable('1TUP')
cmd.hide('all')
cmd.show('sphere', 'name zn')

cmd.show('surface', 'chain A+B+C')
cmd.show('cartoon', 'chain E+F')
cmd.scene('S0', action='store', view=0, frame=0, animate=-1)

cmd.show('cartoon')
cmd.hide('surface')

cmd.scene('S1', action='store', view=0, frame=0, animate=-1)

cmd.hide('cartoon', 'chain A+B+C')
cmd.show('mesh', 'chain A')
cmd.show('sticks', 'chain A+B+C')
cmd.scene('S2', action='store', view=0, frame=0, animate=-1)

cmd.set('ray_trace_mode', 0)
cmd.mset(1, 500)


cmd.frame(0)
cmd.scene('S0')
cmd.mview()
cmd.frame(60)
cmd.set_view((-0.175534308,   -0.331560850,   -0.926960170,
             0.541812420,     0.753615797,   -0.372158051,
             0.821965039,    -0.567564785,    0.047358301,
             0.000000000,     0.000000000, -249.619018555,
             58.625568390,   15.602619171,   77.781631470,
             196.801528931, 302.436492920,  -20.000000000))

cmd.mview()
cmd.frame(90)
cmd.set_view((-0.175534308,   -0.331560850,   -0.926960170,
              0.541812420,    0.753615797,   -0.372158051,
              0.821965039,   -0.567564785,    0.047358301,
              -0.000067875,    0.000017881, -249.615447998,
              54.029174805,   26.956727982,   77.124832153,
             196.801528931,  302.436492920,  -20.000000000))
cmd.mview()
cmd.frame(150)
cmd.set_view((-0.175534308,   -0.331560850,   -0.926960170,
              0.541812420,    0.753615797,   -0.372158051,
              0.821965039,   -0.567564785,    0.047358301,
              -0.000067875,    0.000017881,  -55.406421661,
              54.029174805,   26.956727982,   77.124832153,
              2.592475891,  108.227416992,  -20.000000000))
cmd.mview()
cmd.frame(200)
cmd.scene('S1')
cmd.mview()
cmd.frame(350)
cmd.scene('S1')
cmd.set_view((0.395763457,   -0.173441306,    0.901825786,
              0.915456235,    0.152441502,   -0.372427106,
             -0.072881661,    0.972972929,    0.219108686,
              0.000070953,    0.000013039,  -37.689743042,
             57.748500824,   14.325904846,   77.241867065,
             -15.123448372,   90.511535645,  -20.000000000))

cmd.mview()
cmd.frame(351)
cmd.scene('S2')
cmd.mview()

cmd.frame(500)
cmd.scene('S2')
cmd.mview()
cmd.mplay()
cmd.mpng('p53_1tup')

cmd.quit()
