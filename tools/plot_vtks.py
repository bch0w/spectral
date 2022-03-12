"""
A script to call the pyatoa.visuals.model.Model class and plot standard looking
slices of a model for easy visualization. Contains a set of configurations
to control the look of the figures.

Note:
    One can easily make a .gif to visualize model slices in psuedo-3d by using
    the imagemagick package, with the command. This could be wrapped in
    subprocess but since it's a one liner, it's pretty easy to do quickly.

    $ convert -delay 20 -loop 1 *png model_animated.gif
"""
import os
import sys
import time
from glob import glob
from pyatoa.visuals.vtk_modeler import VTKModeler
from pyatoa.utils.images import imgs_to_pdf


def split_fid(fid):
    """
    VTK file naming needs to be strictly enforced to the template:

    {tag}_{number}_{quantity}.vtk
        *tag: identifier  e.g. 'model', 'gradient', 'update' etc.
        *number: a 4 digit value denoting model number, e.g. 0001
        *quantity: value of the tag, e.g. 'vs', 'vp', 'vs_kernel'
    """
    # Drop directory and suffix
    fid_stripped = os.path.splitext(os.path.basename(fid))[0]
    parts = fid_stripped.split("_")
    tag = parts[0]
    model_number = parts[1]
    # assert int(model_number), f"Model number does not match formatting"
    quantity = " ".join(parts[2:])

    return fid_stripped, tag, model_number, quantity


def save_pdf(fids, fid_out, clean=True):
    """
    Convenience function to collect pngs into a pdf and then remove pngs

    :type fids: list
    :param fids: list of .png files to convert to .pdf
    :type fid_out: str
    :param fid_out: name of output .pdf file
    :type clean: bool
    :pararm clean: delete the original .png files
    """
    # Wait for the image to save, otherwise problems with saving unfinished pngs
    time.sleep(5)
    imgs_to_pdf(fids, fid_out)
    if clean:
        for fid in fids:
            os.remove(fid)


def set_kwargs(vtk_fid, depth_km=None, pct=None, **kwargs):
    """
    Update key word arguments in place based on the model type and depth
    Hard coding for colorscale values etc. should go here. 
    :return: 
    """
    _, tag, model_number, quantity = split_fid(vtk_fid) 

    # Default values are None
    min_max, cbar_title, round_to = None, None, None

    # Determine color kwargs based on tag
    if tag == "gradient":
        cmap = "RdBu"
        reverse = True
        default_range = False
    elif tag in ["log", "update"]:
        cmap = "RdBu"
        reverse = False
        default_range = False
        min_max = [-.25, .25]
        round_to = 1
    elif tag == "poissons":
        cmap = "RdBu"
        reverse = False
        default_range = True
    elif tag in ["divide", "ratio"]:
        cmap = "RdYlBu"
        reverse = True
        round_to = 0.1
        if quantity == "vpvs":
            default_range = False
            min_max = [1.4, 2.25]
        elif quantity == "poissons":
            default_range = False
            min_max = [0.15, 0.4]
    elif tag == "kernel":
        cmap = "RdYlBu"
        reverse = True
        default_range = True
    elif tag == "model":
        cmap = "jet"
        reverse = True
        round_to = 1
        default_range = True
        # if depth_km == "surface":
        #     min_max = {"vp": [2135, 5500], 
        #                "vs": [975, 3250]}[quantity]
        # elif depth_km < 5:
        #     min_max = {"vp": [2135, 5500], 
        #                "vs": [975, 4000]}[quantity]
        # elif 5 <= depth_km < 10:
        #     min_max = {"vp": [3500, 7000], 
        #                "vs": [1750, 4000]}[quantity]
        # elif 10 <= depth_km < 15:
        #     min_max = {"vp": [3500, 7500], 
        #                "vs": [2000, 4000]}[quantity]
        # elif 15 <= depth_km < 30:
        #     min_max = {"vp": [5500, 8000], 
        #                "vs": [2750, 4750]}[quantity]
        # elif 30 <= depth_km < 40:
        #     min_max = {"vp": [6000, 8750], 
        #                "vs": [3000, 5000]}[quantity]
        # elif 40 <= depth_km <= 50:
        #     min_max = {"vp": [7000, 9250], 
        #                "vs": [3500, 5000]}[quantity]
        # else:
        #     default_range = True
    else:
        raise NotImplementedError(tag)

    # Update the kwargs
    kwargs["cmap"] = cmap
    kwargs["reverse"] = reverse
    kwargs["default_range"] = default_range
    kwargs["min_max"] = min_max
    kwargs["round_to"] = round_to
    kwargs["cbar_title"] = f"{tag}\n{quantity}_{model_number}"

    return kwargs


def call_vtk_modeler(vtk_fid, src_fid, rcv_fid, coast_fid, show, path_out,
                     depths, slices, make_pdf=True, **kwargs):
    """
    Call the VTKModeler class to plot models
    """
    fid, _, _, _ = split_fid(vtk_fid) 

    # Make the output directory
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    # Initiate the class with some preset keyword arguments
    vm = VTKModeler(num_clabels=3, num_colors=30, scale_axes=1E-3,
                    zero_origin=True, xlabel="E [km]", ylabel="N [km]",
                    zlabel="Z [km]", figsize=(1000, 1000), offscreen=not show,
                    **kwargs
                    )

    # Load in the VTK files
    vm.load(fid=vtk_fid, src_fid=src_fid, rcv_fid=rcv_fid, coast_fid=coast_fid)

    # Create depth slices and combine into a single PDF
    save_fids = []
    for depth_km in depths:
        save_fid = os.path.join(path_out, f"{fid}_Z_{depth_km}.png")
        if os.path.exists(save_fid):
            print(f"{save_fid} already exists, skipping...")
            continue
        save_fids.append(save_fid)

        vm.kwargs.update(set_kwargs(vtk_fid, depth_km=depth_km, **kwargs))
        vm.depth_slice(depth_km=depth_km, show=show, save=save_fid)

    if make_pdf:
        save_pdf(save_fids, os.path.join(path_out, f"{fid}_Z.pdf"))

    # Create cross-sections for X and Y axes, combine into separate PDFs
    for axis in ["X", "Y"]:
        save_fids = []
        for pct in slices:
            save_fid = os.path.join(path_out,
                                    f"{fid}_{axis}_{pct * 1E2}pct.png")
            if os.path.exists(save_fid):
                print(f"{save_fid} already exists, skipping...")
                continue
            save_fids.append(save_fid)

            vm.kwargs.update(set_kwargs(vtk_fid, pct=pct, **kwargs))
            vm.cross_section(axis=axis, pct=pct, show=show, save=save_fid)

        if save_fids and make_pdf:
            save_pdf(save_fids, os.path.join(path_out, f"{fid}_{axis}.pdf"))

def pdf(fid): 
    """
    MAIN:
    Create PDFs of depth slices and cross sections for each model, gradient etc. 
    Addresses each model individually, deletes original .png files
    """
    # Set slices here
    slices = []  # should be in percentages, e.g. .25
    depths = ["surface", 5, 10, 15]

    # Auxiliary files
    src_fid = "./srcs.vtk"
    rcv_fid = "./rcvs.vtk"
    coast_fid = "/Users/Chow/subduction/data/carto/coastline/coast.npy"

    # Additional parameters
    path_out = f"./figures"
    show = False

    for vtkfid in glob(fid):
        call_vtk_modeler(vtk_fid=vtkfid, src_fid=src_fid, rcv_fid=rcv_fid,
                         coast_fid=coast_fid, path_out=path_out, show=show,
                         slices=slices, depths=depths, make_pdf=False
                         )

def four_banger():
    """
    MAIN:
    Create a tiled image showing four different models, gradients etc. for each
    depth slice. e.g. show Vp, Vs, dVp and dVs all in a single image, useful
    for quick assessments of multiple model parameters

    warning, this is pretty gross coding
    """
    from PIL import Image
    from numpy import array 

    # Set Parameters here
    depths = ["surface", 15]
    slices = []

    # !!! These need to be 'glob'able
    choices = ["model/model_??17_vp*", 
               "model/model_??17_vs*",
               "model/update_??17_vp*",
               "model/update_??17_vs*"
               ]

    # Auxiliary files
    src_fid = "./srcs.vtk"
    rcv_fid = "./rcvs.vtk"
    coast_fid = "/Users/Chow/subduction/data/carto/coastline/coast.npy"

    # We need a temp hold directory to store output file
    scratch = "./scratch"
    path_out = "./figures/panel"

    for fid in [scratch, path_out]:
        if not os.path.exists(fid):
            os.makedirs(fid)

    # Loop on the first list as a 'master list' and derive other fids from that
    # This ensures that all the file names are pre-sorted before making images
    choice_list = [sorted(glob(_)) for _ in choices]
    sorted_choice_list = [[], [], [], []]
    for c1 in choice_list[0]:
        n_mod = os.path.basename(c1).split("_")[1]  # e.g. 0001
        for i in range(len(choice_list)):
            # Find the relevant file based on matching model number
            cx = [_ for _ in choice_list[i] if n_mod in _]
            if not cx:
                print(f"Missing file {n_mod} for {choices[i]}")
                sorted_choice_list[i].append("")
            else:
                sorted_choice_list[i] += cx


    # Loop through the sorted lists for each model number and create the 4 banga
    for sorted_choices in array(sorted_choice_list, dtype=str).T:
        if "" in sorted_choices:
            print(f"Missing files for {choices[0]}, skipping...")
            continue
        # Generate all slices for each choice
        for choice in sorted_choices:
            call_vtk_modeler(vtk_fid=choice, src_fid=src_fid, rcv_fid=rcv_fid,
                             coast_fid=coast_fid, path_out=scratch, show=False,
                             slices=slices, depths=depths, make_pdf=False,
                             )

    # Create four banger using PIL - list sorting should be maintained
    # We need to get the names of the outputted images in order w.r.t each other
    glob_choices = [os.path.join(scratch, _.split("/")[-1]) for _ in choices]
    img_choices = array([sorted(glob(_)) for _ in glob_choices]).T

    for img_list in img_choices:
        # Modified from pyatoa.utils.images.tile_imgs
        images = []
        for fid in img_list:
            images.append(Image.open(fid).convert("RGBA"))
       
        # Widths and heights should match
        widths, heights = zip(*(i.size for i in images))
        im_out = Image.new(mode="RGBA", 
                            size=(max(widths) * 2, max(heights) * 2))

        # Hardcode positions based on images sizes, not very elegant  :/
        # 0 1 
        # 2 3
        im_out.paste(images[0], (0, 0))                 # top left
        im_out.paste(images[1], (images[0].width, 0))   # top right
        im_out.paste(images[2], (0, images[0].height))  # bottom left
        im_out.paste(images[3], (images[2].width, images[2].height))  # btm left

        # Create unique identifier for output file
        _, tag, model_number, quantity = split_fid(img_list[0])
        quantity = "_".join(quantity.split(" ")[1:])

        fid_out = os.path.join(path_out,
                               f"panel_{model_number}_{quantity}.png")
        im_out.save(fid_out)

def single():
    """
    No frills no saving just plot a single slice or cross section for a single
    file and show, useful for quick-looking models
    """
    # Set your arguments here
    vtk_fid = "beta_kernel_smooth.vtk"
    kwargs = {"cmap": "RdBu",
              "reverse": True,
              "default_range": False,
              # "min_max": [-.15, .15],
              # "round_to": None,
              # "cbar_title": "Update Vs_0014"
              }

    # Pick what you want to plot here
    depth_km = 5
    x_pct = None
    y_pct = None

    # Auxiliary files
    src_fid = None # "./srcs.vtk"
    rcv_fid = None # "./rcvs.vtk"
    coast_fid = "/Users/Chow/subduction/data/carto/coastline/coast.npy"

    # Initiate the class with some preset keyword arguments
    vm = VTKModeler(num_clabels=3, num_colors=30, scale_axes=1E-3,
                    zero_origin=True, xlabel="E [km]", ylabel="N [km]",
                    zlabel="Z [km]", figsize=(1000, 1000), offscreen=False,
                    **kwargs
                    )

    # Load in the VTK files
    vm.load(fid=vtk_fid, src_fid=src_fid, rcv_fid=rcv_fid, coast_fid=coast_fid)
    
    if depth_km:
        vm.depth_slice(depth_km=depth_km, show=True, save=False)
    if x_pct:
        vm.cross_section(axis="X", pct=x_pct, show=True, save=False)
    if y_pct:
        vm.cross_section(axis="Y", pct=y_pct, show=True, save=False)

def make_all():
    """
    MAIN:
    Create PDFs of depth slices and cross sections for each model, gradient etc.
    Addresses each model individually, deletes original .png files
    """
    # Set slices here
    slices = []  # should be in percentages, e.g. .25
    depths = ["surface", 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50]

    # File identifier list
    vtk_fids = glob(f"./model/*0011*.vtk") + glob(f"./gradient/*0011*.vtk")

    # Auxiliary files
    src_fid = "./srcs.vtk"
    rcv_fid = "./rcvs.vtk"
    coast_fid = "/Users/Chow/subduction/data/carto/coastline/coast.npy"

    # Additional parameters
    path_out = f"./figures"
    show = True

    # Call the modeler for each fid
    for vtkf in vtk_fids:
        call_vtk_modeler(vtk_fid=vtkf, src_fid=src_fid, rcv_fid=rcv_fid,
                         coast_fid=coast_fid, path_out=path_out, show=show,
                         slices=slices, depths=depths, make_pdf=True
                         )

if __name__ == "__main__":
    single()
    a=1/0
    if len(sys.argv) > 1:
        if sys.argv[1] == "pdf":
            pdf(fid=sys.argv[2])
        elif sys.argv[1] == "one":
            single()
        elif sys.argv[1] == "4banger":
            four_banger()
        else:
            raise NotImplementedError
    else:
        make_all()

