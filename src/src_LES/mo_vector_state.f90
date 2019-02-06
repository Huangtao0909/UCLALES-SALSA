MODULE mo_vector_state
  USE classFieldArray
  USE mo_structured_datatypes  
  IMPLICIT NONE

  SAVE

  TYPE(FloatArray3d), TARGET :: a_up, a_uc, a_ut
  TYPE(FloatArray3d), TARGET :: a_vp, a_vc, a_vt
  TYPE(FloatArray3d), TARGET :: a_wp, a_wc, a_wt
  

  CONTAINS

    SUBROUTINE setVectorVariables(Vector,outputlist,nzp,nxp,nyp)
      TYPE(FieldArray), INTENT(inout) :: Vector
      CHARACTER(len=*), INTENT(in) :: outputlist(:)
      INTEGER, INTENT(in) :: nzp,nxp,nyp
      
      REAL :: zeros3d(nzp,nxp,nyp)
      
      CLASS(*), POINTER :: pipeline_p, pipeline_c, pipeline_t
      
      zeros3d = 0.
      
      a_up = FloatArray3d(zeros3d,store=.TRUE.)
      a_uc = FloatArray3d(zeros3d,store=.TRUE.)
      a_ut = FloatArray3d(zeros3d,store=.TRUE.)
      pipeline_p => a_up
      pipeline_t => a_ut
      pipeline_c => a_uc
      CALL Vector%newField("uwind","Zonal wind vector","m/s",'mttt',         &
                           ANY(outputlist == 'uwind'),in_p_data=pipeline_p,  &
                           in_t_data=pipeline_t,in_c_data=pipeline_c,        &
                           in_group='Vector')

      a_vp = FloatArray3d(zeros3d,store=.TRUE.)
      a_vc = FloatArray3d(zeros3d,store=.TRUE.)      
      a_vt = FloatArray3d(zeros3d,store=.TRUE.)
      pipeline_p => a_vp
      pipeline_t => a_vt
      pipeline_c => a_vc
      CALL Vector%newField("vwind","Meridional wind vector","m/s",'tmtt',         &
                           ANY(outputlist == 'vwind'),in_p_data=pipeline_p,  &
                           in_t_data=pipeline_t,in_c_data=pipeline_c,        &
                           in_group='Vector')
      
      a_wp = FloatArray3d(zeros3d,store=.TRUE.)
      a_wc = FloatArray3d(zeros3d,store=.TRUE.)
      a_wt = FloatArray3d(zeros3d,store=.TRUE.)
      pipeline_p => a_wp
      pipeline_t => a_wt
      pipeline_c => a_wc
      CALL Vector%newField("wwind","Vertical wind vector","m/s",'ttmt',         &
                           ANY(outputlist == 'wwind'),in_p_data=pipeline_p,  &
                           in_t_data=pipeline_t,in_c_data=pipeline_c,        &
                           in_group='Vector')      
                        
    END SUBROUTINE setVectorVariables


    
END MODULE mo_vector_state
