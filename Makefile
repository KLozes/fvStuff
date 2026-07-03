TARGET  = fv_up7
SRC     = fv_up7_ducros_2d.cu

NVCC    = nvcc
ARCH    = native
CFLAGS  = -O3 -arch=$(ARCH) --expt-relaxed-constexpr
LIBS    = -lm

.PHONY: all clean run run_fr run_rt run_keep run_keep_lmv

all: $(TARGET) fr_lob rt_dg dg1 keep_fv

$(TARGET): $(SRC)
	$(NVCC) $(CFLAGS) -o $@ $< $(LIBS)

# Default run: N=400
run: $(TARGET)
	./$(TARGET) 400

FR_SRC  = fr_lobatto_2d.cu

# FR/Lobatto DG targets — set ORDER=1..3 (default 3)
FR_ORDER ?= 3
fr_lob: $(FR_SRC)
	$(NVCC) $(CFLAGS) -DORDER=$(FR_ORDER) -o $@ $< $(LIBS)

run_fr: fr_lob
	./fr_lob 50

RT_SRC = rt_dg_euler_2d.cu

# RT0/P0 mixed DG targets
rt_dg: $(RT_SRC)
	$(NVCC) $(CFLAGS) -o $@ $< $(LIBS)

run_rt: rt_dg
	./rt_dg 400

# Degree-1 Gauss-collocation DG (16 DOF/cell, all DOFs on Gauss, pos. limiter)
# Riemann flux is compile-time: set RIEMANN=RIE_HLLC (default) | RIE_BARMAT | RIE_CENTRAL
RT1_SRC = dg1_euler_2d.cu
dg1: $(RT1_SRC)
	$(NVCC) $(CFLAGS) -o $@ $< $(LIBS)

run_dg1: dg1
	./dg1 200

KEEP_SRC = keep_fv_euler_2d.cu

# Collocated FV with KEEP (ECKEP entropy-conserving + KEP, ES scheme)
keep_fv: $(KEEP_SRC)
	$(NVCC) $(CFLAGS) -o $@ $< $(LIBS)

run_keep: keep_fv
	./keep_fv 200

run_keep_lmv: keep_fv
	./keep_fv 200 lmv 0.1

clean:
	rm -f $(TARGET) fr_lob rt_dg keep_fv sod_*.ppm fr_sod_*.ppm
