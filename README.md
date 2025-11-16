# conductance-measurement

paper name: overcoming the space clamp effect: Reliable recovery of local and effective synaptic conductances of neurons.


我们在四种不同类型的带复杂树突结构和离子通道的神经元模型上说明空间钳效应，并比较了我们的方法和传统方法在恢复局部和等效突触电导的效果。包括L5PC， FS， hippocampus PC和purkinje神经元。其中正文的结果主要是L5PC的结果。由于程序类似，我主要就L5PC的代码为例进行说明：

space clamp effect 说明： space_clamp_error.ipynb (Fig 1A);

deviation between the local synaptic current and measured somatic synaptic current: deviarion_I_syn_local_syn.ipynb (Fig 1B);

recovery of the temporal mean of the local synaptic conductance: local_conductance.ipynb (Fig 1F-II);

recovery of the time constants of the local synaptic conductance: time_constants_est.ipynb (Fig 2); 

Interaction between the injected current and synaptic current: interaction_between_inj_syn.ipynb (Fig 3);

effective conductances estimation:  eff_cond_est.ipynb for a pair of EI inputs (Fig 4C);
                                       
                                   multiple_efff_con.ipynb for multiple E and I inputs (Fig 4D);

paper link: https://www.pnas.org/doi/full/10.1073/pnas.2512294122
