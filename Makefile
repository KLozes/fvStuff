TARGET  = fv_up7
SRC     = fv_up7_ducros_2d.cu

NVCC    = nvcc
ARCH    = native
CFLAGS  = -O3 -arch=$(ARCH) --expt-relaxed-constexpr
LIBS    = -lm

.PHONY: all clean run run_fr run_rt run_p1 run_ned run_sbp run_carb run_keep run_keep_lmv

all: $(TARGET) fr_lob rt_dg rt_dg_dbl sbp_dg p1_dg ned_dg dg1 keep_fv

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

# Double-precision RT0 (low-Mach diagnostics baseline)
rt_dg_dbl: $(RT_SRC)
	$(NVCC) $(CFLAGS) -DUSE_DOUBLE -o $@ $< $(LIBS)

# Strong-form SBP-SAT finite-difference variant: 2-pt Lobatto COLLOCATION of
# the momentum volume flux (instead of exact integration).  With this switch
# the RT0/P0 DG is algebraically identical to du/dt = -Df + M^-1 B (f - f*).
sbp_dg: $(RT_SRC)
	$(NVCC) $(CFLAGS) -DUSE_DOUBLE -DSBP_COLLOC -o $@ $< $(LIBS)

run_sbp: sbp_dg
	./sbp_dg 100 lmv 0.1

P1_SRC = p1_dg_euler_2d.cu

# FULL-LINEAR (P1) momentum / P0 mixed DG — DOUBLE precision.
# Experiment: does full P1 momentum keep low-Mach preservation like RT0?
p1_dg: $(P1_SRC)
	$(NVCC) $(CFLAGS) -o $@ $< $(LIBS)

# Gresho / low-Mach vortex in double precision
run_p1: p1_dg
	./p1_dg 100 lmv 0.1

# NEDELEC (curl-conforming / rotated-RT0) momentum: same source, -DNEDELEC
# projects out the divergence slopes (mxs=mys=0) → div-free, vorticity-carrying.
ned_dg: $(P1_SRC)
	$(NVCC) $(CFLAGS) -DNEDELEC -o $@ $< $(LIBS)

run_ned: ned_dg
	./ned_dg 100 lmv 0.1

CARB_SRC = carb_fv_2d.cu

# Quirk odd-even / carbuncle test — plain 1st-order Godunov FV.
# carb_fv (HLLC, carbuncle-prone) vs carb_fv_hll (HLL, carbuncle-free control)
carb_fv: $(CARB_SRC)
	$(NVCC) $(CFLAGS) -o $@ $< $(LIBS)

carb_fv_hll: $(CARB_SRC)
	$(NVCC) $(CFLAGS) -Driemann_n=hll_n -o $@ $< $(LIBS)

run_carb: carb_fv carb_fv_hll
	./carb_fv 400
	./carb_fv_hll 400

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
	rm -f $(TARGET) fr_lob rt_dg rt_dg_dbl sbp_dg p1_dg ned_dg carb_fv carb_fv_hll dg1 keep_fv sod_*.ppm fr_sod_*.ppm
