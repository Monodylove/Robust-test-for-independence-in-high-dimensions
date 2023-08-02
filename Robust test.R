library(Rcpp)
sourceCpp("Robust_test.cpp")

t1 <- Sys.time()
sN1 <- c(4,8,16,32,64,128,256)
sM1 <- c(4,8,16,32,64,128,256,512)
sN <- c(4,8,16,32)
sM <- c(4,8,16,32,64)
Table1_S_Nm <- RESULT(sN,sM,1)
Table1_S_Nm 
Table1_T_rho <- RESULT(sN,sM,2)
Table1_T_rho 
Table1_T_nm <- RESULT(sN1,sM1,3)
Table1_T_nm 
Table1_T_3 <- RESULT(sN1,sM1,4)
Table1_T_3 
Table4_S_Nm <- RESULT(sN,sM,1,rho = 0.1)
Table4_S_Nm
Table4_T_rho <- RESULT(sN,sM,2,rho = 0.1)
Table4_T_rho
Table4_T_nm <- RESULT(sN1,sM1,3,rho = 0.1)
Table4_T_nm
Table4_T_3 <- RESULT(sN1,sM1,4,rho = 0.1)
Table4_T_3
t2 <- Sys.time()
print(t2-t1)

t3 <- Sys.time()
RESULT(32,128,1,nsim = 500)
RESULT(64,128,1,nsim = 500)
RESULT(128,128,1,nsim = 500)
t4 <- Sys.time()
print(t4-t3)

t5 <- Sys.time()
RESULT(32,128,3,nsim = 500)
RESULT(64,128,3,nsim = 500)
RESULT(128,128,3,nsim = 500)
t6 <- Sys.time()
print(t6-t5)
