import numpy as np


def plot_lcp_nt(plt_n, plt_t, lemke, element_null, nodes, frame_num=None, scale_lcp=None):
    """
    Plots interaction forces and mutual displacements along the normal and tangent to the contact zone
    :param plt_n: matplotlib.pyplot object or axis for plotting contact along the normal to the contact zone
    :param plt_t: matplotlib.pyplot object or axis for plotting contact along the tangent to the contact zone
    :param lemke: results of solving LCP. LCP.lemke.Lemke() class object
    :param element_null: null-elements in scheme. FEM.element_null.ElementNullContainer() class object
    :param nodes: nodes is scheme. FEM.scheme.NodeContainer() class object
    :param frame_num: Used for animation. To get access for lemke.zn_anim[frame_num] or zt_anim[frame_num] ect.
    :param scale_lcp: scale for the mutual displacements (because often interaction forces are bigger)
    :return: None
    """
    if frame_num is None:
        zn, zt, xn, xt = lemke.zn, lemke.zt, lemke.xn, lemke.xt
    else:
        zn, zt, xn, xt = (lemke.zn_anim[frame_num], lemke.zt_anim[frame_num],
                          lemke.xn_anim[frame_num], lemke.xt_anim[frame_num])
    x_range = [nodes[element.EN[0]].x for element in element_null if element.orientation == 'n']
    if scale_lcp is None:
        if np.max(zn) > 0.00000001 and np.max(xn) > 0.0001:  # some smaller value will make plots unreadable
            scale_lcp = int(np.max(xn) / np.max(zn))
        else:
            scale_lcp = 1
    # plot for normal
    plt_n.bar(x_range, xn, color='red', width=0.2)  # plot bars
    for x, y in zip(x_range, xn):  # plot values
        plt_n.text(x, y, str("%.2f" % y), color='lightcoral')
    plt_n.bar(x_range, zn * scale_lcp, color='blue', width=0.2)
    for x, y in zip(x_range, zn):
        plt_n.text(x, y * scale_lcp, str("%.4f" % y), color='cornflowerblue')
    # plot for tangent
    plt_t.bar(x_range, xt, color='red', width=0.2)
    for x, y in zip(x_range, xt):  # plot values
        plt_t.text(x, y, str("%.2f" % y), color='lightcoral')
    plt_t.bar(x_range, zt * scale_lcp, color='blue', width=0.2)
    for x, y in zip(x_range, zt):
        plt_t.text(x, y * scale_lcp, str("%.4f" % y), color='cornflowerblue')




