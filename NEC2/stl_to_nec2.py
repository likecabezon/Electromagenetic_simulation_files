import argparse
import numpy as np
from stl import mesh
import collections
import datetime
import os # For checking file existence

def stl_to_nec2_converter(stl_filepath, output_nec_filepath,
                          wire_radius, wire_resistance_ohms_per_meter,
                          segment_per_edge, comments, ground_plane_mode,
                          input_unit_scale_factor):
    """
    Converts an STL mesh into a NEC-2 wire-grid geometry, ready for
    manual setup of excitations, frequency sweeps, etc.

    Args:
        stl_filepath (str): Path to the input STL file.
        output_nec_filepath (str): Path to save the NEC-2 input file.
        wire_radius (float): Radius of the NEC-2 wires (in meters).
        wire_resistance_ohms_per_meter (float): Ohms per meter for the wire material.
                                                Use 0.0 for perfectly conducting wires.
        segment_per_edge (int): How many NEC-2 segments to create for each STL edge.
                                 If > 1, edges will be discretized.
        comments (list): A list of strings for comments to include in the NEC-2 file.
        ground_plane_mode (int): Ground plane mode (0=Free Space, 1=Infinite Ground Plane with mirroring,
                                 -1=Infinite Ground Plane without mirroring).
        input_unit_scale_factor (float): Multiplier to convert STL input units to meters.
    """

    try:
        # Load the STL mesh using numpy-stl
        your_mesh = mesh.Mesh.from_file(stl_filepath)
    except Exception as e:
        print(f"Error: Could not load STL file '{stl_filepath}'. Please check the path and file integrity.")
        print(f"Details: {e}")
        return

    # Apply the unit scaling to all vertices immediately after loading
    # your_mesh.vectors is a (N, 3, 3) array where N is number of triangles
    # and each triangle has 3 vertices (x,y,z)
    scaled_vectors = your_mesh.vectors * input_unit_scale_factor
    faces = scaled_vectors
    
    print(f"STL loaded: {len(faces)} faces.")
    print(f"Scaling input dimensions by {input_unit_scale_factor} to convert to meters.")

    unique_edges = set() # Set to store unique edges (tuples of sorted vertex tuples)
    nec_segments_data = [] # List to store data for each NEC-2 segment
    current_gw_tag = 1 # Initialize NEC-2 GW tag (wire tag), starts from 1

    # Iterate through each face (triangle) in the STL mesh
    for face in faces:
        # Convert vertex coordinates of the current triangle to tuples for hashing
        v_coords = [tuple(p) for p in face]

        # Process each of the three edges of the current triangle
        for i in range(3):
            p1_coords = v_coords[i]
            p2_coords = v_coords[(i + 1) % 3] # Get the next vertex (wraps around for the third edge)

            # Create a canonical representation of the edge to ensure uniqueness.
            # Sorting the vertex tuples ensures that (v1, v2) and (v2, v1) are treated as the same edge.
            ordered_edge_coords = tuple(sorted((p1_coords, p2_coords)))

            # If this edge has already been processed, skip it to avoid duplicates
            if ordered_edge_coords in unique_edges:
                continue

            unique_edges.add(ordered_edge_coords) # Add the unique edge to the set

            # Convert tuple coordinates back to numpy arrays for numerical calculations
            p1_np = np.array(p1_coords)
            p2_np = np.array(p2_coords)

            # Calculate the vector and length of the current edge
            edge_vector = p2_np - p1_np
            edge_length = np.linalg.norm(edge_vector)

            # Skip zero-length edges to prevent errors
            if edge_length == 0:
                continue

            # Determine the step vector for discretizing the edge into multiple NEC-2 segments
            step_vector = edge_vector / segment_per_edge

            # Create 'segment_per_edge' number of NEC-2 segments for the current STL edge
            for j in range(segment_per_edge):
                seg_start = p1_np + j * step_vector
                seg_end = p1_np + (j + 1) * step_vector

                # Store the data for the current NEC-2 segment
                nec_segments_data.append({
                    'gw_tag': current_gw_tag,
                    'num_segments_in_gw': 1, # This GW card defines a single segment
                    'p1': seg_start,
                    'p2': seg_end,
                    'radius': wire_radius # Wire radius is already expected in meters
                })
                current_gw_tag += 1 # Increment GW tag for the next segment

    print(f"Extracted {len(unique_edges)} unique edges from STL. Generated {len(nec_segments_data)} NEC-2 segments.")
    print(f"Total number of GW cards will be {len(nec_segments_data)}.")


    # Write the NEC-2 input file (.nec)
    with open(output_nec_filepath, 'w') as f:
        # CM cards for comments and header information
        f.write("CM\n")
        f.write(f"CM NEC-2 Wire-Grid Geometry from STL Mesh - Generated on {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("CM This file contains only the geometry (GW cards) and basic loading (LD cards).\n")
        f.write("CM You need to manually add FR (Frequency), EX (Excitation), RP (Radiation Pattern), etc.\n")
        f.write(f"CM Wire Radius: {wire_radius:.6f} meters\n")
        f.write(f"CM Segments per original STL edge: {segment_per_edge}\n")
        
        # Determine and write input unit information
        if input_unit_scale_factor == 1.0:
            f.write("CM Input STL dimensions were treated as meters.\n")
        elif input_unit_scale_factor == 0.01:
            f.write("CM Input STL dimensions were scaled from centimeters to meters (x0.01).\n")
        elif input_unit_scale_factor == 0.001:
            f.write("CM Input STL dimensions were scaled from millimeters to meters (x0.001).\n")
        elif input_unit_scale_factor == 1e-6:
            f.write("CM Input STL dimensions were scaled from micrometers to meters (x1e-6).\n")
        else:
            f.write(f"CM Input STL dimensions were scaled by a factor of {input_unit_scale_factor} to meters.\n")

        if wire_resistance_ohms_per_meter > 0:
            f.write(f"CM Wire Resistance (applied via LD 1): {wire_resistance_ohms_per_meter:.6f} Ohms/meter\n")
        else:
            f.write("CM Wires are considered perfectly conducting (zero resistance).\n")
        
        # Ground plane mode comment
        if ground_plane_mode == 0:
            f.write("CM Ground Plane Mode: Free Space (GE 0)\n")
        elif ground_plane_mode == 1:
            f.write("CM Ground Plane Mode: Infinite Ground Plane with Current Mirroring (GE 1)\n")
            f.write("CM Ensure your geometry is above the ground plane (typically Y=0 or Z=0).\n")
        elif ground_plane_mode == -1:
            f.write("CM Ground Plane Mode: Infinite Ground Plane without Current Mirroring (GE -1)\n")
            f.write("CM You must manually add a GN (Ground Parameters) card to define the ground properties.\n")
            f.write("CM Example (Perfect Ground): GN 1 0 0 0 0 0 0\n")
            f.write("CM Example (Real Ground): GN 2 0 0 0 SIGMA EPSR DEPTH\n")

        # Add any custom comments provided by the user
        if comments:
            f.write("CM Custom Comments:\n")
            for comment_line in comments:
                f.write(f"CM {comment_line}\n")
        
        # CE card marks the end of comments and start of geometry
        f.write("CE\n")

        # GW cards (Geometry Wire)
        # Format: GW ITAG ISEG X1 Y1 Z1 X2 Y2 Z2 RAD
        for seg in nec_segments_data:
            f.write(f"GW {seg['gw_tag']} {seg['num_segments_in_gw']} "
                    f"{seg['p1'][0]:.6f} {seg['p1'][1]:.6f} {seg['p1'][2]:.6f} "
                    f"{seg['p2'][0]:.6f} {seg['p2'][1]:.6f} {seg['p2'][2]:.6f} "
                    f"{seg['radius']:.6f}\n")

        # GE card (Geometry End)
        # GE IP1 IP2 (IP1=0 for free space, IP1=1 for ground plane with mirroring, IP1=-1 for ground without mirroring)
        f.write(f"GE {ground_plane_mode}\n") 

        # LD 1 card (Wire Resistance) if a non-zero resistance is specified
        if wire_resistance_ohms_per_meter > 0:
            f.write("CM\n")
            f.write("CM Applying uniform resistance to all generated wires via individual LD 1 cards.\n")
            f.write("CM ITAG   ISeg1 ISeg2 R_ohmmeter G L C\n")
            for seg in nec_segments_data:
                 f.write(f"LD 1 {seg['gw_tag']} 1 1 {wire_resistance_ohms_per_meter:.6f} 0 0 0\n")
            f.write("CM\n")

        # Placeholder comments for manual user input
        f.write("CE\n")
        f.write("C  *** MANUAL INPUT REQUIRED BELOW THIS LINE ***\n")
        f.write("C  Add your FR (Frequency), EX (Excitation), RP (Radiation Pattern), etc. here.\n")
        f.write("C  Example for Frequency (1 point at 300 MHz): FR 0 0 1 0 300 0\n")
        f.write("C  Example for Voltage Source (on GW 1, segment 1): EX 0 1 1 0 0 0 0\n")
        f.write("C  Example for Radiation Pattern: RP 0 0 1 1000 90 0 0 0 1 0 0\n")
        f.write("CE\n")
        f.write("EN\n") # End of input deck

    print(f"NEC-2 geometry successfully written to '{output_nec_filepath}'.")
    print("\nNext Steps:")
    print(f"1. Open '{output_nec_filepath}' in a text editor.")
    print("2. Add your desired FR (Frequency), EX (Excitation), RP (Radiation Pattern) and other cards.")
    print("3. Ensure the wire radius and resistance are appropriate for your simulation.")
    print("4. If ground plane mode GE -1 was used, remember to add a GN (Ground Parameters) card manually.")
    print("5. Validate the meshing in a NEC-2 viewer (e.g., xnec2c, 4nec2).")

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description="Converts an STL mesh into a NEC-2 wire-grid geometry.",
        formatter_class=argparse.RawTextHelpFormatter # Allows for multi-line help messages
    )
    # Required positional arguments
    parser.add_argument(
        "stl_file",
        type=str,
        help="Path to the input STL file."
    )
    parser.add_argument(
        "wire_radius",
        type=float,
        help="Radius of the NEC-2 wires in meters. This is critical for wire-grid models (e.g., approx. 0.5 * mesh_cell_size)."
    )

    # Optional arguments
    parser.add_argument(
        "-o", "--output_file",
        type=str,
        default="output_nec_geometry.nec",
        help="Name of the output NEC-2 file (default: output_nec_geometry.nec)."
    )
    parser.add_argument(
        "-r", "--resistance",
        type=float,
        default=0.0,
        help="Resistance of the wires in Ohms per meter (default: 0.0 for perfectly conducting wires).\n"
             "Applied via LD 1 cards for each wire segment. For proper surface impedance for wire-grid terminations,\n"
             "consider manually adding LD 5 cards based on material conductivity and frequency."
    )
    parser.add_argument(
        "-s", "--segments_per_edge",
        type=int,
        default=1,
        help="Number of NEC-2 segments to create for each original STL triangle edge (default: 1).\n"
             "Increase this if your STL edges are electrically too long (e.g., > lambda/10)."
    )
    parser.add_argument(
        "-c", "--comment",
        action="append", # Allows the argument to be specified multiple times
        help="Add custom comment lines to the NEC-2 file. Can be used multiple times.\n"
             "Example: -c 'This is a custom comment' -c 'Another line of comment'."
    )
    parser.add_argument(
        "-g", "--ground_plane",
        type=int,
        default=0,
        choices=[0, 1, -1],
        help="Sets the ground plane mode for the GE card (default: 0).\n"
             "  0: Free Space (GE 0)\n"
             "  1: Infinite Ground Plane with Current Mirroring (GE 1)\n"
             " -1: Infinite Ground Plane without Current Mirroring (GE -1). \n"
             "     Requires manual addition of a GN (Ground Parameters) card."
    )

    # Unit conversion arguments - mutually exclusive group
    unit_group = parser.add_mutually_exclusive_group()
    unit_group.add_argument(
        "-u", "--units",
        type=str,
        choices=['m', 'cm', 'mm', 'um'],
        default='m',
        help="Units of the input STL file. Coordinates will be converted to meters.\n"
             "  'm' (meters): No scaling (default).\n"
             "  'cm' (centimeters): Scale by 0.01.\n"
             "  'mm' (millimeters): Scale by 0.001.\n"
             "  'um' (micrometers): Scale by 1e-6."
    )
    unit_group.add_argument(
        "-f", "--scale_factor",
        type=float,
        help="Direct numerical scale factor to convert input STL dimensions to meters.\n"
             "Example: -f 1e-4 for 0.1mm units, or -f 1000 for kilometers to meters.\n"
             "Mutually exclusive with --units."
    )


    args = parser.parse_args()

    # Determine the effective scale factor
    if args.scale_factor is not None:
        effective_scale_factor = args.scale_factor
    else:
        # Map unit strings to scale factors
        unit_map = {
            'm': 1.0,
            'cm': 0.01,
            'mm': 0.001,
            'um': 1e-6
        }
        effective_scale_factor = unit_map[args.units]


    # --- Create a dummy STL file for demonstration if it doesn't exist ---
    dummy_stl_filename = "example_cylinder.stl"
    # Only create dummy if it's the specified input file AND it doesn't exist
    if args.stl_file == dummy_stl_filename and not os.path.exists(dummy_stl_filename):
        print(f"Creating a dummy STL file '{dummy_stl_filename}' for demonstration purposes...")
        try:
            # Simple open cylinder for demo
            from math import sin, cos, pi
            num_segments_cyl = 20 # Number of facets around the circumference
            cyl_radius = 0.05 # 5 cm radius
            cyl_height = 0.1 # 10 cm height

            vertices = []
            faces = []

            # Vertices for the top and bottom circles
            for i in range(num_segments_cyl):
                angle = 2 * pi * i / num_segments_cyl
                x = cyl_radius * cos(angle)
                y = cyl_radius * sin(angle)
                vertices.append([x, y, cyl_height / 2]) # Top circle
                vertices.append([x, y, -cyl_height / 2]) # Bottom circle
            
            # Side faces (forming the cylindrical surface)
            for i in range(num_segments_cyl):
                v_idx_bottom_curr = i * 2
                v_idx_top_curr = i * 2 + 1
                v_idx_bottom_next = ((i + 1) % num_segments_cyl) * 2
                v_idx_top_next = ((i + 1) % num_segments_cyl) * 2 + 1

                # Two triangles forming a quad
                faces.append([v_idx_bottom_curr, v_idx_top_curr, v_idx_top_next])
                faces.append([v_idx_bottom_curr, v_idx_top_next, v_idx_bottom_next])

            cylinder_mesh = mesh.Mesh(np.zeros(len(faces), dtype=mesh.Mesh.dtype))
            for i, f in enumerate(faces):
                for j in range(3):
                    cylinder_mesh.vectors[i][j] = np.array(vertices[f[j]])
            
            cylinder_mesh.save(dummy_stl_filename)
            print(f"Dummy STL '{dummy_stl_filename}' created.")
        except Exception as e:
            print(f"Error creating dummy STL: {e}")
            print("Please provide your own STL file or fix the dummy creation.")
            exit(1) # Exit if dummy file creation fails

    # Execute the conversion with the parsed arguments
    stl_to_nec2_converter(
        stl_filepath=args.stl_file,
        output_nec_filepath=args.output_file,
        wire_radius=args.wire_radius,
        wire_resistance_ohms_per_meter=args.resistance,
        segment_per_edge=args.segments_per_edge,
        comments=args.comment if args.comment else [],
        ground_plane_mode=args.ground_plane,
        input_unit_scale_factor=effective_scale_factor # Pass the calculated scale factor
    )

    # Optional: No automatic cleanup of dummy STL. User can inspect it.
