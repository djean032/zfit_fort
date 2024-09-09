module laser
  use global
  implicit none

  private

  public :: irradiance, power, W

contains
  function irradiance(phi, w_z, rad)
    real(dp), intent(in) :: phi, w_z, rad 
    real(dp) :: irradiance
    irradiance = (2*phi)/(PI_dp*w_z**2)*exp(-2*rad**2/w_z**2)
  end function irradiance

  function power(phi_pk, t, tau)
    real(dp), intent(in) :: phi_pk, t, tau 
    real(dp) :: power, wid
    wid = tau/(2*(ln(2))**0.5)
    power = phi_pk*exp(-(t/wid)**2)
  end function power

  function W(w0, z, M2, wavelength)
    real(dp), intent(in) :: w0, z, M2, wavelength
    real(dp) :: W
    W = w0*sqrt(1 + (z/((PI_dp*w0**2)/(M2*wavelength)))**2)
  end function W
end module laser

