import numba as nb

@nb.njit
def sorted_eigenvalues(vxx, vyy, vzz, vxy, vyz, vzx):
    #assert v.shape == (3,3)
    a, b, c, d, e, f = vxx, vyy, vzz, vxy, vzx, vyz # v[0,0], v[1,1], v[2,2], v[0,1], v[0,2], v[1,2]

    # Analytic eigenvalues solution of the 3x3 input matrix
    tmp1 = -a**2 + a*b + a*c - b**2 + b*c - c**2 - 3*d**2 - 3*e**2 - 3*f**2
    tmp2 = 2*a**3 - 3*a**2*b - 3*a**2*c - 3*a*b**2 + 12*a*b*c - 3*a*c**2 + 9*a*d**2 + 9*a*e**2 - 18*a*f**2 + 2*b**3 - 3*b**2*c - 3*b*c**2 + 9*b*d**2 - 18*b*e**2 + 9*b*f**2 + 2*c**3 - 18*c*d**2 + 9*c*e**2 + 9*c*f**2 + 54*d*e*f
    tmp3 = np.sqrt((4*tmp1**3 + tmp2**2) + 0j)
    tmp4 = (tmp2 + tmp3) ** (1/3)
    tmp5 = 1/3*(a + b + c)
    tmp6 = 1 + 1j*np.sqrt(3)
    tmp7 = 1 - 1j*np.sqrt(3)
    eigv1 = tmp4/(3*2**(1/3)) - (2**(1/3)*tmp1)/(3*tmp4) + tmp5
    eigv2 = (tmp6*tmp1)/(3*2**(2/3)*tmp4) - (tmp7*tmp4)/(6*2**(1/3)) + tmp5
    eigv3 = (tmp7*tmp1)/(3*2**(2/3)*tmp4) - (tmp6*tmp4)/(6*2**(1/3)) + tmp5

    # Assume the values are real ones and remove the FP rounding errors
    eigv1 = np.real(eigv1)
    eigv2 = np.real(eigv2)
    eigv3 = np.real(eigv3)

    # Sort the eigenvalues using a fast sorting network
    eigv1, eigv2 = min(eigv1, eigv2), max(eigv1, eigv2)
    eigv2, eigv3 = min(eigv2, eigv3), max(eigv2, eigv3)
    eigv1, eigv2 = min(eigv1, eigv2), max(eigv1, eigv2)

    return eigv1, eigv2, eigv3