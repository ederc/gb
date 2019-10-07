# TOPSRCDIR allows to compile this testsuite for a clone (and different revision) of this
# repository. If the environment variable topsrcdir is the top source directory of the clone,
# compile this testsuite as:
#
# make TOPSRCDIR="$topsrcdir"
#
# Finally, if you also leave out TOPSRCDIR and just run 'make' then the testsuite
# is compiled against the current source tree.

AUTOMAKE_OPTIONS = foreign
TOPSRCDIR = $(srcdir)/..
TOPBUILDDIR = $(builddir)/..
DEFINES =
# include TOPBUILDIR for m4ri_config.h
AM_CFLAGS = -I$(TOPSRCDIR) -I$(TOPBUILDDIR) \
						${OPENMP_CFLAGS} ${CL_CFLAGS} ${GMP_CFLAGS}

STAGEDIR := $(realpath -s $(TOPBUILDDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) -lgb ${GMP_LIBS} ${CL_LIBS}
LDFLAGS = -no-install

TESTS =												\
	initialize_basis						\
	free_basis									\
	initialize_pairset					\
	free_pairset								\
	initialize_hash_table				\
	free_hash_table							\
	divmask											\
	normalize_row								\
	input_meta_data							\
	xinput_meta_data						\
	julia_io										\
	pairset_update							\
	select_spairs								\
	qq_example									\
	cyclic_4_drl_15_bit					\
	cyclic_4_drl_17_bit					\
	cyclic_4_h_drl_17_bit				\
	cyclic_4_drl_31_bit					\
	cyclic_4_lex_17_bit					\
	cyclic_5_drl_17_bit					\
	cyclic_5_lex_17_bit					\
	cyclic_6_drl_17_bit_la1_t1	\
	cyclic_6_drl_17_bit_la1_t2	\
	cyclic_6_drl_17_bit_la2_t1	\
	cyclic_6_drl_17_bit_la2_t2	\
	cyclic_6_drl_17_bit_la42_t1	\
	cyclic_6_drl_17_bit_la42_t2	\
	cyclic_6_drl_17_bit_la43_t1	\
	cyclic_6_drl_17_bit_la43_t2	\
	cyclic_6_drl_17_bit_la44_t1	\
	cyclic_6_drl_17_bit_la44_t2	\
	katsura_4_lex								\
	katsura_5_lex								\
	cyclic_7										# \
	# katsura_12_drl_17_bit				\
	# katsura_13_drl_17_bit				\
	# cyclic_9										\
	# pearson_9										\
	# katsura_16									\
	# mayr_42											\
	# bayes_148										\
	# ideal_q											\
	# cyclic_6_lex_17_bit

XFAIL_TESTS = 								\
	xinput_meta_data

check_PROGRAMS = $(TESTS) $(XFAIL_TESTS)
