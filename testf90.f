        PROGRAM TimeSeries_Class_Example 
        USE TimeCoordinate_Class 
        USE TimeSeries_Class 
        IMPLICIT NONE 

        REAL, DIMENSION(10), PARAMETER :: values = 
     +        (/ 10.0, 9.0, 7.0, 5.6, 4.3, 1.0, 0.34, 0.01, 0.009, 
     +        0.001 /) 

        TYPE(TimeCoordinate):: time 
        TYPE(TimeSeries)    :: series 
        INTEGER             :: i 

        CALL initialize(time,name="time",unit="days since 000-01-01") 
        CALL initialize(series,name="temperature",time=time, 
     +        description="Time Series of Temperature") 

         DO i=LBound(values,1),UBound(values,1) 
           CALL addNewTime(series,dble(i)) 
           CALL setValue(series,dble(values(i))) 
         END DO 

         CALL print(series) 

         CALL finalize(series) 
         CALL finalize(time) 
       END PROGRAM TimeSeries_Class_Example
