#!/usr/bin/env python

import imageio
import os

import numpy as np
import matplotlib.pyplot as plt

FILE_POINTS = 'data/points.txt'

list_points = []
function_lines = []

def cubic_coefficients(list_points):
    # Calculate cubic spline coefficients for supplied list of points.
    num_points = len(list_points)

    A = np.zeros((num_points, num_points))
    r = np.zeros((num_points, 1))

    dx = lambda i: list_points[i+1][0] - list_points[i][0]
    dy = lambda i: list_points[i+1][1] - list_points[i][1]

    # Populate matrix and right-side, excluding end-conditions.
    for i in range(1,num_points-1):
        A[i, i-1:i+2] = [dx(i-1),2*(dx(i-1)+dx(i)),dx(i)]
        r[i] = 3*(dy(i)/dx(i)-dy(i-1)/dx(i-1))
    # Set up end-conditions in matrix (Natural Cubic Spline).
    A[0][0] = A[-1][-1] = 1;
    # Set up coefficients matrix.
    coeff = np.zeros((num_points, 3))
    # Calculating values for c_1 -> c_n.
    coeff[:,1] = np.linalg.solve(A, r)[:,0]
    # Calculate the other coefficients.
    for i in range(0,num_points-1):
        coeff[i,2] = (coeff[i+1,1] - coeff[i,1])/(3*dx(i))
        coeff[i,0] = dy(i)/dx(i)-dx(i)*(2*coeff[i,1]+coeff[i+1,1])/3
    # Cut the coefficient matrix down to size.
    return coeff[:-1,:]

def cubic_fit_evaluation(list_points, points_per_segment):
    x1 = []
    y1 = []

    x = lambda i: list_points[i][0]
    y = lambda i: list_points[i][1]

    n = len(list_points)
    coeff = cubic_coefficients(list_points)
    for i in range(n-1): # Iterate over each point.
        xs = np.linspace(x(i), x(i+1), points_per_segment)
        dx = xs - x(i)
        ys = coeff[i,2]*dx;
        ys = (ys+coeff[i,1])*dx;
        ys = (ys+coeff[i,0])*dx+y(i);
        # Add to x1 and y1 lists.
        x1 += list(xs)
        y1 += list(ys)
    return (x1, y1)

def tension_coefficients(list_points, tau, y0, yn):

    num_points = len(list_points)

    A = np.zeros((num_points, num_points))
    r = np.zeros((num_points, 1))

    x = lambda i: list_points[i][0]
    y = lambda i: list_points[i][1]

    h = lambda i: x(i+1) - x(i)

    alpha = lambda i: 1/h(i) - tau/np.sinh(tau*h(i))
    beta = lambda i: (tau*np.cosh(tau*h(i)))/np.sinh(tau*h(i)) - 1/h(i)
    gamma = lambda i: (tau*tau*(y(i+1)-y(i)))/h(i)

    # Populate matrix and right side, excluding end-conditions.
    for i in range(1, num_points-1):
        # Select A, row i and elements [i-1,i+1] from that row.
        A[i, i-1:i+2] = [alpha(i-1), beta(i-1)+beta(i), alpha(i)]
        r[i] = gamma(i)-gamma(i-1)

    n = num_points-1
    # Set end-conditions in matrix.
    A[0][0] = -beta(0)
    A[0][1] =  alpha(0)
    A[n][n-1] = alpha(n-1)
    A[n][n] = beta(n-1)
    # Set end-conditions in right hand side.
    r[0] = tau*tau*y0 - gamma(0)
    r[n] = tau*tau*yn - gamma(n-1)

    # Set up coefficient array.
    coeff = np.zeros((num_points, 1))
    # Solve system for coefficients.
    coeff[:,0] = np.linalg.solve(A, r)[:,0]


    return coeff[:,0]

def tension_fit_evaluation(list_points, tau, points_per_segment, y0, yn):

    x1 = []
    y1 = []

    x = lambda i: list_points[i][0]
    y = lambda i: list_points[i][1]
    h = lambda i: x(i+1) - x(i)

    n = len(list_points)

    coeff = tension_coefficients(list_points, tau, y0, yn)
    z = lambda i: coeff[i]

    pow_tau = tau*tau

    for i in range(n-1): # Since there is x_i+1 and z_i+1.
        xs = np.linspace(x(i), x(i+1), points_per_segment)
        dx = xs - x(i)
        dx_plus = x(i+1) - xs

        ## Calculate left side.
        left_result = ((z(i)*np.sinh(tau*dx_plus)) +
                (z(i+1)*np.sinh(tau*dx)))/(pow_tau*np.sinh(tau*h(i)))

        ## Calculate right side.
        right_result = ((y(i)-z(i)/pow_tau)*dx_plus +
        (y(i+1)-z(i+1)/pow_tau)*dx)/h(i)

        # Add values to correct list.
        x1 += list(xs)
        y1 += list(left_result + right_result)

    return (x1, y1)


def on_click(event):
    global list_points
    global function_lines
    if event.inaxes:
        x, y = event.xdata, event.ydata
        list_points.append((x,y))
        plt.plot(x, y, 'bo')
        xs, ys = cubic_fit_evaluation(list_points, 100)
        function_lines.append(plt.plot(xs, ys, 'b')[0])
        ax = plt.gca()
        if (len(function_lines) > 1):
            old_line = function_lines.pop(0)
            ax.lines.remove(old_line)
            del old_line
        plt.gcf().canvas.draw()
    else:
        print("OUTSIDE!")

def main():

    # Import image.
    img = imageio.imread('figs/p2-car-orig.png')
    fig_height, fig_widht, dep = img.shape


    # Reverse the image, so that it starts on 0 y.
    rev_img = [row for row in reversed(img)]
    img = rev_img

    # Calculate width and height for image.
    aspect = fig_widht/fig_height
    h = 2.2 # inches.
    w = h*aspect;


    def display_car():
        fig = plt.figure(figsize = (w, h))
        fig.tight_layout()
        # Display first image and save.
        plt.imshow(img)
        plt.axis([0, fig_widht, 0, fig_height])
        plt.xticks(range(10,fig_widht,20))
        plt.yticks(range(10,fig_height,20))

        ax = fig.axes[0]
        ax.tick_params(labelsize=8)
        plt.xticks(rotation=270)


        plt.title("Car image")
        plt.xlabel('pixel x-coord')
        plt.ylabel('pixel y-coord')
        fig.savefig("figs/p2-car.pdf", bbox_inches='tight')

        LW_GRID = 0.2

        #plt.grid()
        ## Change properties of grid.
        #ax.grid(color='b', linestyle='--', lw=LW_GRID)

        #for x in range(20,fig_widht,20):
        #    plt.plot([x, x], [0, fig_height], 'b--', lw=LW_GRID)

        #for y in range(20,fig_widht,20):
        #    plt.plot([0, fig_widht], [y, y], 'b--', lw=LW_GRID)

        fig.savefig("figs/p2-car-grid.pdf", bbox_inches='tight')
        return fig

    fig = display_car()

    if not os.path.isfile(FILE_POINTS):
        fig.canvas.callbacks.connect('button_press_event', on_click)
        plt.show()
        with open(FILE_POINTS, 'w') as handle:
            txt = '\n'.join(["{},{}".format(x,y) for (x,y) in list_points])
            handle.write(txt)
    else:
        lines = open(FILE_POINTS).readlines()
        read_points = [(float(x),float(y)) for x,y in [line.split(',') for line in lines]]
        xs, ys = cubic_fit_evaluation(read_points, 100)
        plt.plot([x for x,y in read_points], [y for x,y in read_points], 'bo')
        plt.plot(xs, ys, 'b')

    # Always save updated figure.
    plt.title("Car silhouette: cubic spline.")
    fig.savefig("figs/p2-car-grid-cubic.pdf", bbox_inches='tight')


    ## Task 5.
    print("TASK 5")
    fig2, ax2 = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True,
            figsize=(5.2,2.5))

    # add a big axes, hide frame
    fig2.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel("x")
    plt.ylabel("T(x)")

    x_tension = np.linspace(0, 4*np.pi, 10)
    y_tension = [np.cos(x) for x in x_tension]

    points_tension = list(zip(x_tension, y_tension))
    y0 = 0
    yn = 0

    segment_points = 100
    tau = 0.01
    xs, ys = tension_fit_evaluation(points_tension, tau, segment_points, y0, yn)
    ax2[0][0].plot(xs, ys, x_tension, y_tension, 'bo')
    tau = 0.1
    xs, ys = tension_fit_evaluation(points_tension, tau, segment_points, y0, yn)
    ax2[0][1].plot(xs, ys, x_tension, y_tension, 'bo')
    tau = 3.0
    xs, ys = tension_fit_evaluation(points_tension, tau, segment_points, y0, yn)
    ax2[1][0].plot(xs, ys, x_tension, y_tension, 'bo')
    tau = 10.0
    xs, ys = tension_fit_evaluation(points_tension, tau, segment_points, y0, yn)
    ax2[1][1].plot(xs, ys, x_tension, y_tension, 'bo')
    plt.suptitle("Tension splines with different tension values.")
    fig2.savefig("figs/p2-fit-10-points.pdf", bbox_inches='tight')

    fig2_5 = plt.figure(figsize=(5.2,2.7))
    tau = 0.01
    y0 = -5
    yn = -5
    xs, ys = tension_fit_evaluation(points_tension, tau, segment_points, y0, yn)
    plt.plot(x_tension, y_tension, "ko")
    plt.plot(xs, ys, "k-", label="-5")
    y0 = -1
    yn = -1
    xs, ys = tension_fit_evaluation(points_tension, tau, segment_points, y0, yn)
    plt.plot(xs, ys, "k--", label="-1")
    y0 = 3
    yn = 3
    xs, ys = tension_fit_evaluation(points_tension, tau, segment_points, y0, yn)
    plt.plot(xs, ys, "k:", label="3")
    y0 = 10
    yn = 10
    xs, ys = tension_fit_evaluation(points_tension, tau, segment_points, y0, yn)
    plt.plot(xs, ys, "k-.", label="10")
    plt.legend()
    plt.title("Tension splines with varying beginning- and end-conditions.")
    plt.ylabel("T(x)")
    plt.xlabel("x")
    fig2_5.savefig("figs/p2-fit-10-points-y0.pdf", bbox_inches='tight')

    ## Task 6.
    plt.close('all')
    fig3 = display_car()

    file_tension_points = "data/tension_points.txt"
    tension_points = []
    function_lines = []

    def get_points_tension(event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            tension_points.append((x,y))
            plt.plot(x, y, 'bo')
            xs, ys = tension_fit_evaluation(tension_points, tau, 100, y0, yn)
            function_lines.append(plt.plot(xs, ys, 'b')[0])
            ax = plt.gca()
            if (len(function_lines) > 1):
                old_line = function_lines.pop(0)
                ax.lines.remove(old_line)
                del old_line
            plt.gcf().canvas.draw()
        else:
            print("OUTSIDE!")

    tension_params = [
            (0.1, 1, 0),
            (0.1, 0, 0),
            (0.1, 0.2, -0.3),
            (0.8, 0, 0),
        ]

    if not os.path.isfile(file_tension_points):
        fig3.canvas.callbacks.connect('button_press_event', get_points_tension)
        plt.show()
        with open(file_tension_points, 'w') as handle:
            handle.write("\n".join(["{},{}".format(x,y) for x,y in
                tension_points]))
    else:
        with open(file_tension_points) as handle:
            lines = handle.readlines()
            tension_points = [(float(x),float(y)) for x,y in [line.split(',')
                for line in lines]]
        xs = []
        ys = []
        for i, (tau, y0, yn) in enumerate(tension_params):
            t_xs, t_ys = tension_fit_evaluation(tension_points[i*2:i*2+3], tau, 100, y0, yn)
            xs += t_xs
            ys += t_ys
        plt.plot([x for x,y in tension_points], [y for x,y in
            tension_points], 'bo')
        plt.plot(xs, ys, 'b')

    plt.title("Car silhouette: tension spline.")
    fig3.savefig("figs/p2-car-tension-splines.pdf", bbox_inches='tight')

    file_tension_points = "data/tension_points_extra.txt"
    tension_points = []
    function_lines = []

    tension_params = [
            (1.0, 0, 0),
            (0.01, -0.8, 1),
            (5.0, 0.0, -0.0),
            (0.01, -0.3, 1.5),
            (5.0, 0, 0),
        ]

    if not os.path.isfile(file_tension_points):
        fig3.canvas.callbacks.connect('button_press_event', get_points_tension)
        plt.show()
        with open(file_tension_points, 'w') as handle:
            handle.write("\n".join(["{},{}".format(x,y) for x,y in
                tension_points]))
    else:
        with open(file_tension_points) as handle:
            lines = handle.readlines()
            tension_points = [(float(x),float(y)) for x,y in [line.split(',')
                for line in lines]]
        xs = []
        ys = []
        for i, (tau, y0, yn) in enumerate(tension_params):
            t_xs, t_ys = tension_fit_evaluation(tension_points[i*2:i*2+3], tau, 100, y0, yn)
            xs += t_xs
            ys += t_ys
        plt.plot([x for x,y in tension_points], [y for x,y in
            tension_points], 'bo')
        plt.plot(xs, ys, 'b')

    file_tension_points = "data/tension_points_extra_window.txt"
    tension_points = []
    function_lines = []

    tension_params = [
            (0.001, 1, -0.2),
            (0.01, 0.2, -0.6),
#            (5.0, 0.0, -0.0),
#            (0.01, -0.3, 1.5),
#            (5.0, 0, 0),
        ]

    if not os.path.isfile(file_tension_points):
        fig3.canvas.callbacks.connect('button_press_event', get_points_tension)
        plt.show()
        with open(file_tension_points, 'w') as handle:
            handle.write("\n".join(["{},{}".format(x,y) for x,y in
                tension_points]))
    else:
        with open(file_tension_points) as handle:
            lines = handle.readlines()
            tension_points = [(float(x),float(y)) for x,y in [line.split(',')
                for line in lines]]
        xs = []
        ys = []
        for i, (tau, y0, yn) in enumerate(tension_params):
            t_xs, t_ys = tension_fit_evaluation(tension_points[i*2:i*2+3], tau, 100, y0, yn)
            xs += t_xs
            ys += t_ys
        plt.plot([x for x,y in tension_points], [y for x,y in
            tension_points], 'bo')
        plt.plot(xs, ys, 'b')

    file_tension_points = "data/tension_points_extra_window_lower.txt"
    tension_points = []
    function_lines = []

    tension_params = [
            (0.001, -0.1, 0.2),
#            (0.01, 0.3, -0.6),
#            (5.0, 0.0, -0.0),
#            (0.01, -0.3, 1.5),
#            (5.0, 0, 0),
        ]

    if not os.path.isfile(file_tension_points):
        fig3.canvas.callbacks.connect('button_press_event', get_points_tension)
        plt.show()
        with open(file_tension_points, 'w') as handle:
            handle.write("\n".join(["{},{}".format(x,y) for x,y in
                tension_points]))
    else:
        with open(file_tension_points) as handle:
            lines = handle.readlines()
            tension_points = [(float(x),float(y)) for x,y in [line.split(',')
                for line in lines]]
        xs = []
        ys = []
        for i, (tau, y0, yn) in enumerate(tension_params):
            t_xs, t_ys = tension_fit_evaluation(tension_points[i*2:i*2+3], tau, 100, y0, yn)
            xs += t_xs
            ys += t_ys
        plt.plot([x for x,y in tension_points], [y for x,y in
            tension_points], 'bo')
        plt.plot(xs, ys, 'b')

    plt.title("Car silhouette: tension spline.")
    fig3.savefig("figs/p2-car-tension-splines-extra.pdf", bbox_inches='tight')
    fig3.set_size_inches(10,10)
    fig3.savefig("figs/p2-car-tension-splines-extra-large.pdf",
            bbox_inches='tight')

#
#    fig3 = plt.figure()
#    xs = np.linspace(0,2,100)
#    f = lambda x: (x*np.cosh(x))/np.sinh(x)
#    f2 = lambda x: np.cosh(x)/np.sinh(x)
#    print("beta with tau zero: {}".format(f(0)))
#    plt.plot(xs, [f2(x) for x in xs], label="coth")
#    plt.plot(xs, [np.sinh(x) for x in xs], label="sinh")
#    plt.plot(xs, [np.cosh(x) for x in xs], label="cosh")
#    plt.plot(xs, [f(x) for x in xs], label="beta")
#    plt.legend()
#    plt.show()


#    print(list_points)

#    plt.show()
#
#    x = np.arange(0, 400, 0.1)
#    y = np.sin(x)*250/2
#    y += 250/2
#
#    plt.plot(x, y)
#
#    plt.show()
#
#    plt.figure(figsize=(

if __name__ == "__main__":
    main()

