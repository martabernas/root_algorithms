

def main():
    x0_bis, neval_bis = bisection(opt_fun, [-5, 1], 1e-6)
    x0_sec, neval_sec = secant_root(opt_fun, [-5, 1], 1e-6)
    x0_nwt, neval_nwt = newton_root(opt_fun, opt_fun_deriv, -5, 1e-6)

    print(f"Root of a tested function determined with bisection is: {x0_bis}.\n"
          f"Number of function evaluation: {neval_bis}")
    print(f"Root of a tested function determined with secant method is: {x0_sec}.\n"
          f"Number of function evaluation: {neval_sec}")
    print(f"Root of a tested function determined with Newton method is: {x0_nwt}.\n"
          f"Number of function evaluation: {neval_nwt}")


def bisection(fun, range_list, y_error):
    """Returns the root of given function determined with bisection algorithm and number of function
    evaluations. Takes in function handle (fun), list with searching interval range (range_list)
    and desired precision (y_error). Searching interval must contain sign change. """

    xa = range_list[0]
    xb = range_list[1]
    neval = 0
    xmid = xa
    if (fun(xa) * fun(xb)) < 0:
        while abs(fun(xmid)) > y_error:
            neval = neval + 1
            xmid = (xa + xb) / 2
            fmid = fun(xmid)
            if fmid * fun(xa) < 0:
                xb = xmid
            else:
                xa = xmid
    else:
        print("Search aborted. Interval has to contain a sign change.")

    return xmid, neval


def secant_root(fun, initial_range, y_error):
    """Returns the root of given function determined with secant method and number of function
    evaluations. Takes in function handle (fun), list with initial range (initial_range)
    and desired precision (y_error)."""

    x1 = initial_range[0]
    x2 = initial_range[1]
    x3 = x1
    neval = 0

    while abs(fun(x2)) > y_error:
        neval = neval + 1
        x3 = x2 - fun(x2) * (x2 - x1) / (fun(x2) - fun(x1))
        x1 = x2
        x2 = x3

    return x3, neval


def newton_root(fun, fun_deriv, x_init, y_error):
    """Returns the root of given function determined with newton method and number of function
    evaluations. Takes in handles to function (fun) and its derivative (fun_deriv), starting point
    (x_init) and desired precision (y_error)."""

    x1 = x_init
    x2 = x1
    neval = 0

    while abs(fun(x2)) > y_error:
        neval = neval + 1
        x2 = x1 - fun(x1) / fun_deriv(x1)
        x1 = x2

    return x2, neval


def opt_fun(x):
    return (x - 1) ** 2 - 2


def opt_fun_deriv(x):
    return 2 * (x - 1)


if __name__ == '__main__':
    main()
