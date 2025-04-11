import taskblaster as tb
# from ase.build import nanoribbon


# This is a taskblaster workflow that generates a series of workflows for
# generating armchair graphene nanoribbons (AGNR) of different widths.
@tb.dynamical_workflow_generator_task
def generate_wfs_from_list(input):
    for i in range(1, input+1):
        wf = SubWorkflow(width=i)
        name = f'AGNR_{i}'
        yield name, wf


@tb.workflow
class SubWorkflow:
    width = tb.var()

    @tb.task
    def generate_ribbon(self):
        return tb.node('print_task', number=self.width)


def print_task(number):
    print(f'AGNR_{number} generated')
    return number
