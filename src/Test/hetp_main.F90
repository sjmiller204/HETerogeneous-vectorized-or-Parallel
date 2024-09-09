!######### HETP: AEROSOL THERMODYNAMIC EQUILIBRIUM OF THE ###################
!#########       NH4-SO4-NO3-Na-Cl-Ca-Mg-K SYSTEM         ###################
!
!Copyright (C) 2023  Stefan Miller, Environment and Climate Change Canada
!                    Contact: Stefan.Miller (at) ec.gc.ca
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!############################################################################
! ## HETP Code: SAMPLE MAIN SUBROUTINE TO INTERACT WITH HETP 
!
! ## Copyright 2023, Environment and Climate Change Canada (ECCC)
! ## Written by Stefan Miller
!############################################################################
program hv_main
   use hetp_mod
   implicit none 
!
   integer, parameter   :: dp = kind(0.0d0) 
   integer, parameter   :: sp = kind(0.0e0) 

!  INPUT VARIABLES (total gas + aerosol, input as mol/m3 air)
   real(dp) :: TS          ! Total available sulfate
   real(dp) :: TA          ! Total avaiable ammonium
   real(dp) :: TN          ! Total available nitrate
   real(dp) :: TNa         ! Total available sodium
   real(dp) :: TCl         ! Total available chloride 
   real(dp) :: TCa         ! Total available calcium
   real(dp) :: TK          ! Total available potassium
   real(dp) :: TMg         ! Total available magnesium
   real(dp) :: rh          ! Relative humidity (0-1 scale)
   real(dp) :: temp        ! Air temperature (K)
!
!  OUTPUT VARIABLES (output as mol/m3 air)
   real(dp) :: so4         ! SO4--     (aq) 
   real(dp) :: hso4        ! HSO4-     (aq)
   real(dp) :: caso4       ! CaSO4     (s)
   real(dp) :: nh4         ! NH4+      (aq)
   real(dp) :: nh3         ! NH3       (g)
   real(dp) :: no3         ! NO3-      (aq)
   real(dp) :: hno3        ! HNO3      (g)
   real(dp) :: cl          ! Cl-       (aq)
   real(dp) :: hcl         ! HCl       (g)
   real(dp) :: na          ! Na+       (aq)
   real(dp) :: ca          ! Ca2+      (aq)
   real(dp) :: k           ! K+        (aq)
   real(dp) :: mg          ! Mg2+      (aq)
   real(dp) :: h           ! H+
   real(dp) :: oh          ! OH-
   real(dp) :: lwc         ! Aerosol liquid water
   real(dp) :: frso4       ! Free SO4 
   real(dp) :: frna        ! Free Na
   real(dp) :: frca        ! Free Ca
   real(dp) :: frk         ! Free K
   real(dp) :: frmg        ! Free Mg
   real(dp) :: case_number ! Chemical subspace case number 
!
!  Initalize variables 
   so4    = 0.0_dp
   hso4   = 0.0_dp
   caso4  = 0.0_dp
   nh4    = 0.0_dp
   nh3    = 0.0_dp
   no3    = 0.0_dp
   hno3   = 0.0_dp
   cl     = 0.0_dp
   hcl    = 0.0_dp
   na     = 0.0_dp
   ca     = 0.0_dp
   k      = 0.0_dp
   mg     = 0.0_dp
   h      = 0.0_dp 
   oh     = 0.0_dp
   lwc    = 0.0_dp
   case_number = 0.0_dp
!
!  ## Call HETP for a single set of initial conditions ###################
!  #######################################################################
!  ### Set input data ### umol/m^3 air -> mol/m^3 air ####################
   TS   = 6.122e-2_dp * 1.0e-6_dp             
   TA   = 5.882e-2_dp * 1.0e-6_dp             
   TN   = 3.175e-2_dp * 1.0e-6_dp              
   TNa  = 4.348e-2_dp * 1.0e-6_dp               
   TCl  = 4.110e-2_dp * 1.0e-6_dp             
   TCa  = 1.247e-2_dp * 1.0e-6_dp              
   TK   = 1.279e-2_dp * 1.0e-6_dp              
   TMg  = 4.115e-3_dp * 1.0e-6_dp           
   rh   = 0.90_dp             
   temp = 298.0_dp          
!  #######################################################################
!
!  ### Call HETP
   call mach_hetp_main_15cases(TS, TA, TN, TNa, TCl, TCa, TK, TMg, temp, rh,     &
                               so4, hso4, caso4, nh4, nh3, no3, hno3, cl, hcl,   &
                               na, ca, k, mg, h, oh, lwc, frna, frca, frk, frmg, &
                               frso4, case_number)
!
!  ### Save HETP output to file as columns; column meaning below ###
!  ### Units = mol/m^3 air #########################################
!  *ALWC = aerosol liquid water content 
! [1]  SO4--(aq)   [2]  HSO4-(aq)   [3]  CaSO4(s)   [4]  Free SO4
! [5]  NH4+(aq)    [6]  NH3(g)      [7]  NO3-(aq)   [8]  HNO3(g)
! [9]  Cl-(aq)     [10] HCl(g)      [11] Na+(aq)    [12] Free Na
! [13] Ca++(aq)    [14] Free Ca     [15] K+(aq)     [16] Free K
! [17] Mg++(aq)    [18] Free Mg     [19] ALWC       [20] H+
! [21] OH-         [22] Case Number
   open(unit=51, file='./alloutput_hv.txt', status='replace')
   write(51, 22), so4, hso4, caso4, frso4, nh4, nh3, no3,       &
                  hno3, cl, hcl, na, frna, ca, frca, k,         &
                  frk, mg, frmg, lwc, h, oh, case_number
  22         format(1x, 22(1x, e20.14))

!
!  ### Option to write inputs and outputs to terminal
   !write(*, *), 'HETP inputs configured in HETP/src/Test/hetp_main.F90'
   !write(*, *), '   TS   [mol/m3] = ', TS
   !write(*, *), '   TA   [mol/m3] = ', TA
   !write(*, *), '   TN   [mol/m3] = ', TN
   !write(*, *), '   TNa  [mol/m3] = ', TNa
   !write(*, *), '   TCl  [mol/m3] = ', TCl
   !write(*, *), '   TCa  [mol/m3] = ', TCa
   !write(*, *), '   TK   [mol/m3] = ', TK
   !write(*, *), '   TMg  [mol/m3] = ', TMg
   !write(*, *), '   rh   [frac]   = ', rh
   !write(*, *), '   temp [K]      = ', temp
   !write(*, *), ' '
   !write(*, *), 'HETP outputs'
   !write(*, *), '   so4   : ', so4
   !write(*, *), '   hso4  : ', hso4
   !write(*, *), '   caso4 : ', caso4
   !write(*, *), '   frso4 : ', frso4
   !write(*, *), '   nh4   : ', nh4
   !write(*, *), '   nh3   : ', nh3
   !write(*, *), '   no3   : ', no3
   !write(*, *), '   hno3  : ', hno3
   !write(*, *), '   cl    : ', cl
   !write(*, *), '   hcl   : ', hcl
   !write(*, *), '   na    : ', na
   !write(*, *), '   frna  : ', frna
   !write(*, *), '   ca    : ', ca
   !write(*, *), '   frca  : ', frca
   !write(*, *), '   k     : ', k
   !write(*, *), '   frk   : ', frk
   !write(*, *), '   mg    : ', mg
   !write(*, *), '   frmg  : ', frmg
   !write(*, *), '   lwc   : ', lwc
   !write(*, *), '   h     : ', h
   !write(*, *), '   oh    : ', oh
   !write(*, *), '   case_number : ', case_number
!
stop    
end
