module laser
  use stdlib_constants
  implicit none

  private

  public :: irradiance, power, W

contains
  subroutine irradiance(phi, w_z, r, I)
    real(8), intent(in) :: phi, w_z, r 
    real(8) :: I
    I = (2*phi)/(PI_dp*w_z**2)*exp(-2*r**2/w_z**2)
  end subroutine irradiance

  subroutine power(phi_pk, t, tau, phi)
    real(8), intent(in) :: phi_pk, t, tau
    real(8) :: phi
    phi = phi_pk*exp((-t/tau)**2)
  end subroutine power

  function W(w0, z, M2, wavelength)
    real(8), intent(in) :: w0, z, M2, wavelength
    real(8) :: W
    W = w0*sqrt(1 + (z/((PI_dp*w0**2)/(M2*wavelength)))**2)
  end function W
end module laser

