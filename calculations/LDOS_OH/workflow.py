import taskblaster as tb


@tb.workflow
class Workflow:
    n = tb.var()
    m = tb.var()

    @tb.task
    def generate_ribbon_task(self):
        return tb.node('generate_ribbon',
                       repeat=5,
                       n=self.n,
                       m=self.m)

    @tb.task
    def tag_atoms_task(self):
        return tb.node('tag_atoms',
                       atoms_path=self.generate_ribbon_task,
                       m=self.m)

    @tb.task
    def gs_ribbon_task(self):
        gs_params = {"func": "PBE", "basis": "dzp"}
        return tb.node('gs_calculation',
                       atoms_path=self.tag_atoms_task,
                       params=gs_params,
                       file_name='GS_ribbon.gpw')

    # Edge state
    @tb.task
    def ldos_edge_state_task(self):
        return tb.node('LDOS_func',
                       atoms_path=self.gs_ribbon_task,
                       tag_number=10)

    @tb.task
    def plot_ldos_edge_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_edge_state_task['energies'],
                       ldos=self.ldos_edge_state_task['ldos'],
                       title='LDOS at the interface',
                       filename='pre_relax_interface.png')

    # Middle of the clean part
    @tb.task
    def ldos_clean_state_task(self):
        return tb.node('LDOS_func',
                       atoms_path=self.gs_ribbon_task,
                       tag_number=20)

    @tb.task
    def plot_ldos_clean_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_clean_state_task['energies'],
                       ldos=self.ldos_clean_state_task['ldos'],
                       title='LDOS in the prestine region',
                       filename='pre_relax_prestine.png')

    # Middle of the serated part
    @tb.task
    def ldos_serated_state_task(self):
        return tb.node('LDOS_func',
                       atoms_path=self.gs_ribbon_task,
                       tag_number=30)

    @tb.task
    def plot_ldos_serated_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_serated_state_task['energies'],
                       ldos=self.ldos_serated_state_task['ldos'],
                       title='LDOS in the serated region',
                       filename='pre_relax_serated.png')

    @tb.task
    def add_OH_task(self):
        return tb.node('insert_OH',
                       atoms_path=self.tag_atoms_task)

    @tb.task
    def gs_OH_task(self):
        gs_params = {"func": "PBE", "basis": "dzp"}
        return tb.node('gs_calculation',
                       atoms_path=self.add_OH_task,
                       params=gs_params,
                       file_name='GS_OH.gpaw')

    # Edge state
    @tb.task
    def ldos_OH_pre_relax_task(self):
        return tb.node('LDOS_func',
                       atoms_path=self.gs_OH_task,
                       tag_number=10)

    @tb.task
    def plot_OH_pre_relax_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_OH_pre_relax_task['energies'],
                       ldos=self.ldos_OH_pre_relax_task['ldos'],
                       title='LDOS at the interface with OH pre relaxation',
                       filename='pre_relax_OH_interface.png')

    @tb.task
    def relax_task(self):
        relax_params = {"func": "PBE",
                        "basis": ['szp(dzp)', 'szp(dzp)', 'dzp']}
        return tb.node('multi_relaxation',
                       atoms_path=self.add_OH_task,
                       filename='relaxation.txt',
                       params=relax_params,
                       k_arr=[1, 3, 6])

    # Edge state
    @tb.task
    def ldos_OH_post_relax_task(self):
        return tb.node('LDOS_func',
                       atoms_path=self.relax_task,
                       tag_number=10)

    @tb.task
    def plot_ldos_OH_post_relax_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_OH_post_relax_task['energies'],
                       ldos=self.ldos_OH_post_relax_task['ldos'],
                       title='LDOS at the interface with OH post relaxation',
                       filename='post_relax_OH_interface.png')


def workflow(runner):
    runner.run_workflow(Workflow(n=1, m=1))
