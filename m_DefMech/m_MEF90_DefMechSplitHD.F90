#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechSplitHD,MEF90_DIM)D
#define MEF90_DEFMECHSPLITHD_CONSTRUCTOR MEF90_APPEND(m_MEF90_DefMechSplitHD_Constructor,MEF90_DIM)D
#include "finclude/petscdef.h"

   Use m_MEF90
   Use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
   implicit none

   Type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHSPLITHD
      PetscReal                                        :: gamma
   Contains
      Procedure, pass(self)                            :: EED   => EEDHD
      Procedure, pass(self)                            :: DEED  => DEEDHD
      Procedure, pass(self)                            :: D2EED => D2EEDHD
   End Type

   interface MEF90_DEFMECHSPLITHD
      module procedure MEF90_DEFMECHSPLITHD_CONSTRUCTOR
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DEFMECHSPLITHD_CONSTRUCTOR"
!!!
!!!  
!!!  MEF90_DEFMECHSPLITHD_CONSTRUCTOR: the default constructor for a MEF90_DEFMECHSPLITHD
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Type(MEF90_DEFMECHSPLITHD) Function MEF90_DEFMECHSPLITHD_CONSTRUCTOR(gamma)
      PetscReal,Intent(IN)                             :: gamma

      MEF90_DEFMECHSPLITHD_CONSTRUCTOR%gamma       = gamma
      MEF90_DEFMECHSPLITHD_CONSTRUCTOR%damageOrder = 3
      MEF90_DEFMECHSPLITHD_CONSTRUCTOR%strainOrder = 2
      MEF90_DEFMECHSPLITHD_CONSTRUCTOR%type        = 'MEF90DefMech_unilateralContactTypeHD'
   End Function MEF90_DEFMECHSPLITHD_CONSTRUCTOR



#undef __FUNCT__
#define __FUNCT__ "EEDHD"
!!!
!!!  
!!!  EEDHD: Compute the positive and negative part of the elastic energy density associated with a strain tensor Epsilon
!!!              following the expression from [Amor et al., 2009] W^+ = [sigma(Epsilon)]^+ . Epsilon^+ = sigma(Epsilon) . Epsilon^+
!!!              by orthogality of the masonry projection
!!! [Amor et al., 2009] Amor, H., Marigo, J.-J., and Maurini, C. (2009). Regularized formulation of the variational brittle fracture with unilateral contact: 
!!!                     Numerical experiments. J. Mech. Phys. Solids, 57(8):1209 – 1229.
!!!    
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!    
   Subroutine EEDHD(self,Strain,HookesLaw,EEDPlus,EEDMinus)
      Class(MEF90_DEFMECHSPLITHD),Intent(IN)           :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MATS3D)                                     :: Strain3D
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      PetscReal, Intent(OUT)                           :: EEDPlus,EEDMinus

      PetscErrorCode                                   :: ierr
      Character(len=MEF90_MXSTRLEN)                    :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
      	  ! Only valid for (at least) cubic symmetry
          EEDMinus  = (MEF90_DefMechSplit_SmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%fullTensorLocal%XXXX + 2.0_Kr*HookesLaw%fullTensorLocal%XXYY) / 3.0_Kr) * 0.5_Kr
#if MEF90_DIM==2
         If (HookesLaw%isPlaneStress) Then
            Write(IOBuffer,*) "Hydrostatic-Deviatoric projection not implemented for plane stress anisotropic Hooke's Law: "//__FUNCT__//"\n"
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr)
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer,ierr)
         Else
            Strain3D%XX = Strain%XX
            Strain3D%YY = Strain%YY
            Strain3D%ZZ = 0.0_Kr
            Strain3D%XY = Strain%XY
            Strain3D%XZ = 0.0_Kr
            Strain3D%YZ = 0.0_Kr
            EEDPlus     = ((MEF90_DefMechSplit_SmoothPositiveSquare( trace(Strain),self%gamma) * (HookesLaw%fullTensorLocal%XXXX + 2.0_Kr*HookesLaw%fullTensorLocal%XXYY) / 3.0_Kr) + &
                        & ((HookesLaw%fullTensor3D * (Strain3D - (trace(Strain)/3.0_Kr)*MEF90MatS3DIdentity)) .dotP. (Strain3D - (trace(Strain)/3.0_Kr)*MEF90MatS3DIdentity)) ) * 0.5_Kr
                     ! ((HookesLaw%fullTensor3D  * Strain3D) .dotP. Strain3D) * 0.5_Kr - EEDMinus
         End If
#elif MEF90_DIM==3
     	  EEDPlus   = ((HookesLaw * Strain) .dotP. Strain) * 0.5_Kr - EEDMinus
#endif
!#if MEF90_DIM==2
!         EEDMinus  = (MEF90_DefMechSplit_SmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%fullTensor%XXXX + HookesLaw%fullTensor%YYYY + &
!                     & 2.*HookesLaw%fullTensor%XXYY) / MEF90_DIM**2)  * 0.5_Kr
!#elif MEF90_DIM==3
!         EEDMinus  = (MEF90_DefMechSplit_SmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%fullTensor%XXXX + HookesLaw%fullTensor%YYYY + HookesLaw%fullTensor%ZZZZ + &
!                     & 2.*HookesLaw%fullTensor%XXYY + 2.*HookesLaw%fullTensor%XXZZ + 2.*HookesLaw%fullTensor%YYZZ) / MEF90_DIM**2)  * 0.5_Kr
!#endif
         ! Write(IOBuffer,*) "Hydrostatic-Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         ! Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr)
         ! SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer,ierr)
      !End If
      Else 
         EEDMinus  = MEF90_DefMechSplit_SmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM)  * 0.5_Kr
         EEDPlus   = ((HookesLaw * Strain) .dotP. Strain) * 0.5_Kr - EEDMinus
      End If
      !print *,"K = ", (HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM)
      !print *,"k = ", (HookesLaw%fullTensorLocal%XXXX + 2.0_Kr*HookesLaw%fullTensorLocal%XXYY) / 3.0_Kr
      !print *, "Trace = ", MEF90_DefMechSplit_SmoothPositiveSquare(-trace(Strain),self%gamma)
      !print *,"Strain = ", Strain, "  Stress = ", HookesLaw * Strain
      !print *,"EEDPLUS = ", EEDPLUS, "  EEDMINUS = ", EEDMINUS
      !print *,"trace+ = ", (MEF90_DefMechSplit_SmoothPositiveSquare( trace(Strain),self%gamma) * (HookesLaw%fullTensorLocal%XXXX + 2.0_Kr*HookesLaw%fullTensorLocal%XXYY) / 3.0_Kr)
      !print *,"ETOT = ", ((HookesLaw * Strain) .dotP. Strain) * 0.5_Kr
   End Subroutine

#undef __FUNCT__
#define __FUNCT__ "DEEDHD"
!!!
!!!  
!!!  DEEDHD: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine DEEDHD(self,Strain,HookesLaw,DEEDPlus,DEEDMinus)
      Class(MEF90_DEFMECHSPLITHD),Intent(IN)           :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MATS3D)                                     :: Strain3D
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      Type(MEF90_MATS),Intent(OUT)                     :: DEEDPlus,DEEDMinus
      Type(MATS3D)                                     :: DEEDPlus3D

      PetscErrorCode                                   :: ierr
      Character(len=MEF90_MXSTRLEN)                    :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         ! Only valid for (at least) cubic symmetry
         DEEDMinus  = ((MEF90_DefMechSplit_DSmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%fullTensorLocal%XXXX + 2.0_Kr*HookesLaw%fullTensorLocal%XXYY) / 3.0_Kr) * 0.5_Kr) * MEF90_MATS_IDENTITY
#if MEF90_DIM==2
         If (HookesLaw%isPlaneStress) Then
            Write(IOBuffer,*) "Hydrostatic-Deviatoric projection not implemented for plane stress anisotropic Hooke's Law: "//__FUNCT__//"\n"
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr)
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer,ierr)
         Else
            Strain3D%XX = Strain%XX
            Strain3D%YY = Strain%YY
            Strain3D%ZZ = 0.0_Kr
            Strain3D%XY = Strain%XY
            Strain3D%XZ = 0.0_Kr
            Strain3D%YZ = 0.0_Kr
            DEEDPlus3D  = ((MEF90_DefMechSplit_DSmoothPositiveSquare( trace(Strain),self%gamma) * (HookesLaw%fullTensorLocal%XXXX + 2.0_Kr*HookesLaw%fullTensorLocal%XXYY) / 3.0_Kr) * 0.5_Kr * MEF90MatS3DIdentity) + &
                      & (HookesLaw%fullTensor3D * (Strain3D - (trace(Strain)/3.0_Kr)*MEF90MatS3DIdentity))
            DEEDPlus%XX = DEEDPlus3D%XX
            DEEDPlus%YY = DEEDPlus3D%YY
            DEEDPlus%XY = DEEDPlus3D%XY
         End If
#elif MEF90_DIM==3
     	  DEEDPlus   = (HookesLaw * Strain) - DEEDMinus
#endif
!#if MEF90_DIM==2
!         DEEDMinus = ((-MEF90_DefMechSplit_DSmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%fullTensor%XXXX + HookesLaw%fullTensor%YYYY + &
!                     & 2.*HookesLaw%fullTensor%XXYY) / MEF90_DIM**2 ) * 0.5_Kr) * MEF90_MATS_IDENTITY
!#elif MEF90_DIM==3
!         DEEDMinus = ((-MEF90_DefMechSplit_DSmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%fullTensor%XXXX + HookesLaw%fullTensor%YYYY + HookesLaw%fullTensor%ZZZZ + &
!                     & 2.*HookesLaw%fullTensor%XXYY + 2.*HookesLaw%fullTensor%XXZZ + 2.*HookesLaw%fullTensor%YYZZ) / MEF90_DIM**2 ) * 0.5_Kr) * MEF90_MATS_IDENTITY
!#endif
         !Write(IOBuffer,*) "Hydrostatic-Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         !Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr)
         !SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer,ierr)
      Else
         DEEDMinus = (-MEF90_DefMechSplit_DSmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM) * 0.5_Kr) * MEF90_MATS_IDENTITY
         DEEDPlus  =  (HookesLaw * Strain) - DEEDMinus
      End If
      Call PetscLogFlops(6._pflop,ierr);CHKERRQ(ierr)
   End Subroutine

#undef __FUNCT__
#define __FUNCT__ "D2EEDHD"
!!!
!!!  
!!!  D2EEDHD: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine D2EEDHD(self,Strain,HookesLaw,D2EEDPlus,D2EEDMinus)
      Class(MEF90_DEFMECHSPLITHD),Intent(IN)           :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      Type(MEF90_HOOKESLAW),Intent(OUT)                :: D2EEDPlus,D2EEDMinus
      Type(MEF90HookesLaw3D)                           :: D2EEDPlus3D

      PetscErrorCode                                   :: ierr
      Character(len=MEF90_MXSTRLEN)                    :: IOBuffer
      PetscReal                                        :: ddplus, ddminus

      D2EEDPlus%fullTensor  = 0.0_Kr
      D2EEDMinus%fullTensor = 0.0_Kr
      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
#if MEF90_DIM==2
         If (HookesLaw%isPlaneStress) Then
            Write(IOBuffer,*) "Hydrostatic-Deviatoric projection not implemented for plane stress anisotropic Hooke's Law: "//__FUNCT__//"\n"
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr)
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer,ierr)
         Else
            ddminus = (MEF90_DefMechSplit_D2SmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%fullTensorLocal%XXXX + 2.0_Kr*HookesLaw%fullTensorLocal%XXYY) / 3.0_Kr) * 0.5_Kr
            D2EEDMinus%fullTensor%XXXX = ddminus
            D2EEDMinus%fullTensor%YYYY = ddminus
            D2EEDMinus%fullTensor%XXYY = ddminus
            ddplus = (MEF90_DefMechSplit_D2SmoothPositiveSquare( trace(Strain),self%gamma) * (HookesLaw%fullTensorLocal%XXXX + 2.0_Kr*HookesLaw%fullTensorLocal%XXYY) / 3.0_Kr) * 0.5_Kr
            D2EEDPlus3D%fullTensor%XXXX = ddplus
            D2EEDPlus3D%fullTensor%YYYY = ddplus
            D2EEDPlus3D%fullTensor%ZZZZ = ddplus
            D2EEDPlus3D%fullTensor%XXYY = ddplus
            D2EEDPlus3D%fullTensor%XXZZ = ddplus
            D2EEDPlus3D%fullTensor%YYZZ = ddplus
            D2EEDPlus3D%fullTensor = D2EEDPlus3D%fullTensor + Tens4OS3DXTens4OS3D(HookesLaw%fullTensor3D, (MEF90Tens4OS3DIdentity - ((1.0_Kr/3.0_Kr) * (MEF90MatS3DIdentity .oDot. MEF90MatS3DIdentity))))
            D2EEDPlus%fullTensor%XXXX = D2EEDPlus3D%fullTensor%XXXX
            D2EEDPlus%fullTensor%YYYY = D2EEDPlus3D%fullTensor%YYYY
            D2EEDPlus%fullTensor%XYXY = D2EEDPlus3D%fullTensor%XYXY
            D2EEDPlus%fullTensor%XXYY = D2EEDPlus3D%fullTensor%XXYY
            D2EEDPlus%fullTensor%YYXY = D2EEDPlus3D%fullTensor%YYXY
            D2EEDPlus%fullTensor%XXXY = D2EEDPlus3D%fullTensor%XXXY
         End If
#elif MEF90_DIM==3
         ddminus = (MEF90_DefMechSplit_D2SmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%fullTensor%XXXX + HookesLaw%fullTensor%YYYY + HookesLaw%fullTensor%ZZZZ + &
              & 2.*HookesLaw%fullTensor%XXYY + 2.*HookesLaw%fullTensor%XXZZ + 2.*HookesLaw%fullTensor%YYZZ) / MEF90_DIM**2) * 0.5_Kr
         D2EEDMinus%fullTensor%XXXX = ddminus
         D2EEDMinus%fullTensor%YYYY = ddminus
         D2EEDMinus%fullTensor%ZZZZ = ddminus
         D2EEDMinus%fullTensor%XXYY = ddminus
         D2EEDMinus%fullTensor%XXZZ = ddminus
         D2EEDMinus%fullTensor%YYZZ = ddminus
         D2EEDPlus = HookesLaw - D2EEDMinus
#endif
!#if MEF90_DIM==2
!         dd = (MEF90_DefMechSplit_D2SmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%fullTensor%XXXX + HookesLaw%fullTensor%YYYY + &
!              & 2.*HookesLaw%fullTensor%XXYY) / MEF90_DIM**2) * 0.5_Kr
!         D2EEDMinus%fullTensor%XXXX = dd
!         D2EEDMinus%fullTensor%YYYY = dd
!         D2EEDMinus%fullTensor%XXYY = dd
!#elif MEF90_DIM==3
!         dd = (MEF90_DefMechSplit_D2SmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%fullTensor%XXXX + HookesLaw%fullTensor%YYYY + HookesLaw%fullTensor%ZZZZ + &
!              & 2.*HookesLaw%fullTensor%XXYY + 2.*HookesLaw%fullTensor%XXZZ + 2.*HookesLaw%fullTensor%YYZZ) / MEF90_DIM**2) * 0.5_Kr
!         D2EEDMinus%fullTensor%XXXX = dd
!         D2EEDMinus%fullTensor%YYYY = dd
!         D2EEDMinus%fullTensor%ZZZZ = dd
!         D2EEDMinus%fullTensor%XXYY = dd
!         D2EEDMinus%fullTensor%XXZZ = dd
!         D2EEDMinus%fullTensor%YYZZ = dd      
!#endif
      Else         
         !Write(IOBuffer,*) "Hydrostatic-Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         !Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr)
         !SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer,ierr)
      !End If

         !D2EEDPlus%fullTensor  = 0.0_Kr
         !D2EEDMinus%fullTensor = 0.0_Kr
         D2EEDPlus%Type    = MEF90HookesLawTypeIsotropic
         D2EEDMinus%Type   = MEF90HookesLawTypeIsotropic
         D2EEDMinus%lambda = MEF90_DefMechSplit_D2SmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM) * 0.5_Kr
         D2EEDMinus%mu     = 0.0_Kr
#if MEF90_DIM == 2
         D2EEDMinus%isPlaneStress = HookesLaw%isPlaneStress
         D2EEDMinus%YoungsModulus = 2.0_Kr * D2EEDMinus%mu * (1.0_Kr + D2EEDMinus%PoissonRatio)
         D2EEDMinus%BulkModulus   = D2EEDMinus%lambda + D2EEDMinus%mu
         If (HookesLaw%isPlaneStress) Then
            D2EEDMinus%PoissonRatio  = D2EEDMinus%lambda / (D2EEDMinus%lambda + D2EEDMinus%mu) * 0.5_Kr
            Call PetscLogFlops(14._pflop,ierr);CHKERRQ(ierr)
         Else
            D2EEDMinus%PoissonRatio  = D2EEDMinus%lambda / (D2EEDMinus%lambda + 2.0_Kr * D2EEDMinus%mu) * 0.5_Kr
            Call PetscLogFlops(17._pflop,ierr);CHKERRQ(ierr)
         End If
#else
         D2EEDMinus%PoissonRatio  = D2EEDMinus%lambda / (D2EEDMinus%lambda + D2EEDMinus%mu) * 0.5_Kr
         D2EEDMinus%YoungsModulus = D2EEDMinus%mu * (3.0_Kr * D2EEDMinus%lambda + 2.0_Kr * D2EEDMinus%mu) / (D2EEDMinus%lambda + D2EEDMinus%mu)
         D2EEDMinus%BulkModulus   = D2EEDMinus%lambda + D2EEDMinus%mu * 2.0_Kr / 3.0_Kr
         Call PetscLogFlops(19._pflop,ierr);CHKERRQ(ierr)
#endif
         D2EEDPlus = HookesLaw - D2EEDMinus
      End If
   End Subroutine
End Module MEF90_APPEND(m_MEF90_DefMechSplitHD,MEF90_DIM)D
