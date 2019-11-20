C $Header:  $
C $Name:  $

C--   COMMON /PCEL_MIX/ Integer valued parameters used by the model.
C     interViscAr_pCell :: account for partial-cell in interior vert. visc.
C     interDiffKr_pCell :: account for partial-cell in vertical diffusivity
C     pCellMix_select   :: select option to enhance mixing near surface & bottom
C                          unit digit: near bottom ; tens digit: near surface
C                          with digit =0 : disable ;
C                         = 1 : increases mixing linearly with recip_hFac
C                         = 2,3,4 : increases mixing by recip_hFac^(2,3,4)
C     pCellMix_maxFac  :: maximum enhanced mixing factor
C     pCellMix_delR    :: thickness criteria   for too thin partial-cell
C     pCellMix_diffKr  :: vertical diffusivity for too thin partial-cell
C     pCellMix_viscAr  :: vertical viscosity   for too thin partial-cell

      COMMON /PCELL_MIX_L/
     &        interDiffKr_pCell, interViscAr_pCell
      LOGICAL interDiffKr_pCell
      LOGICAL interViscAr_pCell

      COMMON /PCELL_MIX_I/
     &        pCellMix_select
      INTEGER pCellMix_select

      COMMON /PCELL_MIX_R/
     &        pCellMix_maxFac, pCellMix_delR,
     &        pCellMix_diffKr, pCellMix_viscAr
      _RL pCellMix_maxFac
      _RL pCellMix_delR
      _RL pCellMix_diffKr(Nr)
      _RL pCellMix_viscAr(Nr)

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
