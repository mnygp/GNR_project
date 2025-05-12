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

    # Edge state
    @tb.task
    def ldos_edge_state_task(self):
        return tb.node('gs_and_LDOS',
                       atoms_path=self.tag_atoms_task,
                       tag_number=10)

    @tb.task
    def plot_ldos_edge_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_edge_state_task['energies'],
                       ldos=self.ldos_edge_state_task['ldos'],
                       title='LDOS at the interface',
                       filename='interface.png')

    # Middle of the clean part
    @tb.task
    def ldos_clean_state_task(self):
        return tb.node('gs_and_LDOS',
                       atoms_path=self.tag_atoms_task,
                       tag_number=20)

    @tb.task
    def plot_ldos_clean_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_clean_state_task['energies'],
                       ldos=self.ldos_clean_state_task['ldos'],
                       title='LDOS in the pristine region',
                       filename='prestine.png')

    # Middle of the serated part
    @tb.task
    def ldos_serated_state_task(self):
        return tb.node('gs_and_LDOS',
                       atoms_path=self.tag_atoms_task,
                       tag_number=30)

    @tb.task
    def plot_ldos_serated_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_serated_state_task['energies'],
                       ldos=self.ldos_serated_state_task['ldos'],
                       title='LDOS in the serrated region',
                       filename='serrated.png')

    @tb.task
    def add_OH_task(self):
        return tb.node('insert_OH',
                       atoms_path=self.tag_atoms_task)

    # Edge state
    @tb.task
    def ldos_OH_task(self):
        return tb.node('gs_and_LDOS',
                       atoms_path=self.add_OH_task,
                       tag_number=10)

    @tb.task
    def plot_OH_task(self):
        return tb.node('plot_ldos',
                       energies=self.ldos_OH_task['energies'],
                       ldos=self.ldos_OH_task['ldos'],
                       title='LDOS at the interface with OH',
                       filename='OH_interface.png')


def workflow(runner):
    runner.run_workflow(Workflow(n=1, m=1))
