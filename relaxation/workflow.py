import taskblaster as tb


@tb.workflow
class Workflow:

    @tb.task
    def generate_ribbon_task(self):
        return tb.node('ribbon_path', N=5, identifier='S', n=1, m=1)

    @tb.task
    def relax_ribbon_PW_task(self):
        params = {
            "func": "PBE",
            "PW_cut": 600
        }
        return tb.node('single_PW_relax', atoms_path=self.generate_ribbon_task,
                       filename='calc_out_PW.txt',
                       params=params)

    @tb.task
    def relax_ribbon_dzp_task(self):
        params = {
            "func": "PBE",
            "basis": "dzp"
        }
        return tb.node('single_LCAO_relax', filename='calc_out_LCAO.txt',
                       atoms_path=self.generate_ribbon_task, params=params)

    @tb.task
    def relax_ribbon_LCAO_multi_step_task(self):
        params = {
            "func": "PBE",
            "basis_list": ["szp", "szp", "dzp"]
        }
        k_arr = [3, 6, 6]
        return tb.node('multi_relaxation', filename='calc_out_multi.txt',
                       atoms_path=self.generate_ribbon_task,
                       params=params,
                       k_arr=k_arr)


def workflow(runner):
    return runner.run_workflow(Workflow())
