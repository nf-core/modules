from nf_core.components.update import ComponentUpdate


class SubworkflowUpdate(ComponentUpdate):
    def __init__(
        self,
        pipeline_dir,
        force=False,
        prompt=False,
        sha=None,
        update_all=False,
        show_diff=None,
        save_diff_fn=None,
        update_deps=False,
        remote_url=None,
        branch=None,
        no_pull=False,
        limit_output=False,
    ):
        super().__init__(
            pipeline_dir,
            "subworkflows",
            force,
            prompt,
            sha,
            update_all,
            show_diff,
            save_diff_fn,
            update_deps,
            remote_url,
            branch,
            no_pull,
            limit_output,
        )
