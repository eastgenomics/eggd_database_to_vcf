"""
Functions to communicate with DNAnexus
"""

import os
import dxpy


def download_input_file(remote_file) -> str:
    """
    Download given input file with same name as file in project
    Function from vcf_qc.py from eggd_vcf_qc

    Parameters
    ----------
    remote_file : str
        Name or ID of input CSV file in DNAnexus

    Returns
    -------
    str
        name of locally downloaded file
    """
    local_name = dxpy.describe(remote_file).get("name")
    dxpy.bindings.dxfile_functions.download_dxfile(
        dxid=remote_file, filename=local_name
    )

    return local_name


def upload_output_file(outfile) -> dict:
    """
    Upload output file to set folder in current project
    Function from vcf_qc.py from eggd_vcf_qc

    Parameters
    ----------
    outfile : str
        name of file to upload
    """
    output_project = os.environ.get("DX_PROJECT_CONTEXT_ID")
    output_folder = (
        dxpy.bindings.dxjob.DXJob(os.environ.get("DX_JOB_ID"))
        .describe()
        .get("folder", "/")
    )
    print(f"\nUploading {outfile} to {output_project}:{output_folder}")

    dxpy.set_workspace_id(output_project)
    dxpy.api.project_new_folder(
        output_project, input_params={"folder": output_folder, "parents": True}
    )

    url_file = dxpy.upload_local_file(
        filename=outfile,
        folder=output_folder,
        wait_on_close=True,
    )

    return dxpy.dxlink(url_file)
