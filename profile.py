import subprocess
import re
import sys

def run_test(n_proc):
    print(f"Running with {n_proc} processes...")
    cmd = ["mpirun", "-np", str(n_proc), "./microPF"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        return result.stdout
    except subprocess.TimeoutExpired:
        print(f"Timeout with {n_proc} processes")
        return ""
    except Exception as e:
        print(f"Error with {n_proc} processes: {e}")
        return ""

def parse_timer_output(output):
    # Look for the timer table
    lines = output.split('\n')
    table_start = -1
    for i, line in enumerate(lines):
        if "+---------------------------------------------+------------+------------+" in line:
            table_start = i
            break
    
    if table_start == -1:
        return {}

    data = {}
    # Skip header lines
    for line in lines[table_start+3:]:
        if "+---------------------------------------------+------------+------------+" in line:
            break
        parts = line.split('|')
        if len(parts) > 3:
            section = parts[1].strip()
            time_str = parts[2].strip()
            try:
                # Format is usually "1.234s" or similar, deal.II formatting varies
                # But TimerOutput::summary usually gives:
                # | Section | time | % |
                # Let's assume standard deal.II output
                if time_str.endswith('s'):
                    time_val = float(time_str[:-1])
                else:
                    time_val = float(time_str.split()[0]) # Handle "1.23 s"
                data[section] = time_val
            except:
                pass
    return data

def main():
    # Build first
    subprocess.run(["make", "-j8"], check=True)

    results = {}
    procs = [1, 2, 4, 8]
    
    for n in procs:
        output = run_test(n)
        data = parse_timer_output(output)
        results[n] = data
        
        # Also get total wall time from the output if possible, or sum of sections
        total_time = 0
        if "Total wallclock time elapsed since start" in output:
             # This string might be different depending on deal.II version
             pass

    print("\nScaling Results (Time in seconds):")
    print(f"{'Section':<30} | {'1 Proc':<10} | {'2 Proc':<10} | {'4 Proc':<10} | {'8 Proc':<10}")
    print("-" * 80)
    
    sections = set()
    for n in procs:
        sections.update(results[n].keys())
    
    for section in sorted(sections):
        row = f"{section:<30}"
        for n in procs:
            val = results[n].get(section, "N/A")
            if isinstance(val, float):
                row += f" | {val:<10.4f}"
            else:
                row += f" | {val:<10}"
        print(row)

if __name__ == "__main__":
    main()
