import math
import numpy as np

def generate_cantenna_nec_file(filename, radius, length, max_wire_length, can_metal_thickness, freq, monopole_pos_resolution, pos, len, monopole_rad):
    """
    Generate an NEC file for simulating a cantenna.

    Parameters:
    filename : str
        The name of the output NEC file.
    radius : float
        The radius of the cylindrical cantenna base in meters.
    length : float
        The height or length of the cylindrical cantenna in meters.
    max_wire_length : float
        The maximum length of each wire segment in meters.
    can_metal_thickness : float
        The thickness of the metal used for the cantenna in meters.
    freq : float
        The operating frequency of the cantenna in Hertz.
    monopole_pos_resolution : float
        The resolution for positioning the monopole along the height of the cantenna.
    pos : float
        The desired position of the monopole along the height of the cantenna in meters.
    len : float
        The length of the monopole in meters.
    monopole_rad : float
        The radius of the monopole in meters.

    Returns:
    None
    """

    num_points = int(np.ceil(2 * np.pi * radius / max_wire_length))
    angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
    
    num_circumferences_base = int(np.ceil(radius/max_wire_length))
    base_circumference_radius = np.linspace(0, radius, num_circumferences_base)

    num_heights = int(np.ceil(length/max_wire_length))
    inter_height_num = int(np.ceil(max_wire_length/monopole_pos_resolution))
    heights = np.linspace(0, length, num_heights * inter_height_num + 1)

    aprox_pos = heights[(np.abs(heights - pos)).argmin()]

    
    with open(filename, 'w') as nec_file:
        nec_file.write("CM Cantenna\n")
        nec_file.write("CE\n")

        wire_number = 1

        # Generate wires for the base
        for i in range(1, num_circumferences_base):
            for j in range(num_points):
                x1 = base_circumference_radius[i] * math.cos(angles[j])
                y1 = base_circumference_radius[i] * math.sin(angles[j])
                x2 = base_circumference_radius[i-1] * math.cos(angles[j])
                y2 = base_circumference_radius[i-1] * math.sin(angles[j])
                nec_file.write(f"GW {wire_number} 1 {x1:.9f} {y1:.9f} 0 {x2:.9f} {y2:.9f} 0 {can_metal_thickness:.9f}\n")
                wire_number += 1

                x1 = base_circumference_radius[i] * math.cos(angles[j])
                y1 = base_circumference_radius[i] * math.sin(angles[j])
                x2 = base_circumference_radius[i] * math.cos(angles[j-1])
                y2 = base_circumference_radius[i] * math.sin(angles[j-1])
                nec_file.write(f"GW {wire_number} 1 {x1:.9f} {y1:.9f} 0 {x2:.9f} {y2:.9f} 0 {can_metal_thickness:.9f}\n")
                wire_number += 1


        # Generate wires for the body
        for i in range(1, num_heights + 1):
            for j in range(num_points):
                if j != 0:
                    
                    x1 = radius * math.cos(angles[j])
                    y1 = radius * math.sin(angles[j])
                    z1 = heights[(i) * inter_height_num]
                    x2 = radius * math.cos(angles[j])
                    y2 = radius * math.sin(angles[j])
                    z2 = heights[(i-1) * inter_height_num]
                    nec_file.write(f"GW {wire_number} 1 {x1:.9f} {y1:.9f} {z1:.9f} {x2:.9f} {y2:.9f} {z2:.9f} {can_metal_thickness:.9f}\n")
                    wire_number += 1

                x1 = radius * math.cos(angles[j])
                y1 = radius * math.sin(angles[j])
                z1 = heights[i * inter_height_num]
                x2 = radius * math.cos(angles[j-1])
                y2 = radius * math.sin(angles[j-1])
                z2 = heights[i * inter_height_num]
                nec_file.write(f"GW {wire_number} 1 {x1:.9f} {y1:.9} {z1:.9f} {x2:.9f} {y2:.9f} {z2:.9f} {can_metal_thickness:.9f}\n")
                wire_number += 1


        for i in range(num_heights * inter_height_num):
            x1 = radius * math.cos(angles[0])
            y1 = radius * math.sin(angles[0])
            z1 = heights[i+1]
            x2 = radius * math.cos(angles[0])
            y2 = radius * math.sin(angles[0])
            z2 = heights[i]
            nec_file.write(f"GW {wire_number} 1 {x1:.9f} {y1:.9} {z1:.9f} {x2:.9f} {y2:.9f} {z2:.9f} {can_metal_thickness:.9f}\n")
            wire_number += 1

        nec_file.write(f"GW {wire_number} 14 {radius:.9f} 0 {aprox_pos:.9f} {(radius-len):.9f} 0 {aprox_pos:.9f} {monopole_rad:.9f}\n")
        nec_file.write(f"GE 0\n")
        
        nec_file.write(f"GN	-1\n")

        nec_file.write(f"EK\n")

        nec_file.write(f"EX 0 {wire_number} 1 0 1 0 0 0\n")

        nec_file.write(f"FR 0 0 0 0 {freq/1e6:.3f} 0\n")
        
        nec_file.write(f"XQ\n")

        nec_file.write(f"EN")

    print(f"NEC file '{filename}' generated successfully.")

# Example usage
generate_cantenna_nec_file('cantenna_demo.nec', 0.05, 0.2, 0.01, 1e-4, 2450e6, 0.005, 0.042, 0.027, 5e-4)
