import numpy as np
import mayavi.mlab as mlab
import time



def plot_u(u, x, y, t, n, do = False, save_plot = False):
    """User action function for plotting."""
    if do == False:
        return 0
    else:
        if t[n] == 0:
            time.sleep(2)
    # Mayavi visualization
        mlab.clf()
        extent1 = (0, 20, 0, 20,-2, 2)
        s = mlab.surf(x , y, u,
                          colormap='Blues',
                          extent=extent1)
        mlab.axes(s, color=(.7, .7, .7), extent=extent1,
                      ranges=(0, 10, 0, 10, -1, 1),
                      xlabel='', ylabel='', zlabel='',
                      x_axis_visibility=False,
                      z_axis_visibility=False)
        mlab.outline(s, color=(0.7, .7, .7), extent=extent1)
        mlab.text(6, -2.5, '', z=-4, width=0.14)
        mlab.colorbar(object=None, title=None,
                          orientation='horizontal',
                          nb_labels=None, nb_colors=None,
                          label_fmt=None)
        mlab.title('Waves t=%g' % t[n])
        mlab.view(142, -72, 50)
        f = mlab.gcf()
        camera = f.scene.camera
        camera.yaw(0)

        time.sleep(0) # pause between frames
        if save_plot:
            filename = 'tmp_%04d.png' % n
            mlab.savefig('./animations/'+filename)  # time consuming!
        mlab.show()
