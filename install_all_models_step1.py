import os

installationDirectory = os.path.abspath(os.path.dirname(__file__))

filelist = ['connectome/control/B1border/86_CDK_20041214_BAC_run5_soma_Hay2013_B1border_apic_rec_scaled_diameters.param', \
            'connectome/control/B2border/86_CDK_20041214_BAC_run5_soma_Hay2013_B2border_apic_rec_scaled_diameters.param', \
            'connectome/control/B3border/86_CDK_20041214_BAC_run5_soma_Hay2013_B3border_apic_rec_scaled_diameters.param', \
            'connectome/control/C1border/86_CDK_20041214_BAC_run5_soma_Hay2013_C1border_apic_rec_scaled_diameters.param', \
            'connectome/control/C2center/86_CDK_20041214_BAC_run5_soma_Hay2013_C2center_apic_rec_scaled_diameters.param', \
            'connectome/control/C3border/86_CDK_20041214_BAC_run5_soma_Hay2013_C3border_apic_rec_scaled_diameters.param', \
            'connectome/control/D1border/86_CDK_20041214_BAC_run5_soma_Hay2013_D1border_apic_rec_scaled_diameters.param', \
            'connectome/control/D2border/86_CDK_20041214_BAC_run5_soma_Hay2013_D2border_apic_rec_scaled_diameters.param', \
            'connectome/control/D3border/86_CDK_20041214_BAC_run5_soma_Hay2013_D3border_apic_rec_scaled_diameters.param', \
            'evoked_activity/reduced_model/C2_evoked_UpState_EXC_generic_INH_generic_C2center.param', \
            'celltype_PSTH/ongoing_activity_celltype_template_exc_conductances_fitted.param']

for path in filelist:
    path = os.path.join(installationDirectory, path)
    with open(path, 'r') as in_:
        newStr = in_.read().replace('[INSTALLATION_DIR]', installationDirectory)
    with open(path, 'w') as out_:
        out_.write(newStr)
