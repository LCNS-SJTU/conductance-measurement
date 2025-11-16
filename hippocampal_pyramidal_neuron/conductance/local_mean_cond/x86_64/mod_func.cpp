#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern "C" void _h_reg(void);
extern "C" void _ICaL_reg(void);
extern "C" void _ICaT_reg(void);
extern "C" void _kadist_reg(void);
extern "C" void _kaprox_reg(void);
extern "C" void _kdrca1_reg(void);
extern "C" void _na3_reg(void);
extern "C" void _nax_reg(void);
extern "C" void _synexp_reg(void);

extern "C" void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/hippocampal_pyramidal_neuron/mod/h.mod\"");
    fprintf(stderr, " \"/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/hippocampal_pyramidal_neuron/mod/ICaL.mod\"");
    fprintf(stderr, " \"/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/hippocampal_pyramidal_neuron/mod/ICaT.mod\"");
    fprintf(stderr, " \"/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/hippocampal_pyramidal_neuron/mod/kadist.mod\"");
    fprintf(stderr, " \"/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/hippocampal_pyramidal_neuron/mod/kaprox.mod\"");
    fprintf(stderr, " \"/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/hippocampal_pyramidal_neuron/mod/kdrca1.mod\"");
    fprintf(stderr, " \"/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/hippocampal_pyramidal_neuron/mod/na3.mod\"");
    fprintf(stderr, " \"/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/hippocampal_pyramidal_neuron/mod/nax.mod\"");
    fprintf(stderr, " \"/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/hippocampal_pyramidal_neuron/mod/synexp.mod\"");
    fprintf(stderr, "\n");
  }
  _h_reg();
  _ICaL_reg();
  _ICaT_reg();
  _kadist_reg();
  _kaprox_reg();
  _kdrca1_reg();
  _na3_reg();
  _nax_reg();
  _synexp_reg();
}
