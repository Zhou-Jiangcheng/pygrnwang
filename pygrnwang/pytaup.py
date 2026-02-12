import os
import sys
import platform

import numpy as np
import jpype.imports  # use jpype to call java class

if platform.system() == "Windows":
    jar_path = os.path.join(sys.exec_prefix, 'Scripts', 'TauP.jar')
else:
    jar_path = os.path.join(sys.exec_prefix, 'bin', 'TauP.jar')

if not os.path.exists(jar_path):
    print(f"TauP.jar not found in {sys.exec_prefix}/bin or {sys.exec_prefix}/Scripts")
    print("Please install the TauP toolkit and ensure TauP.jar is in the correct directory.")
    sys.exit(1)

if not jpype.isJVMStarted():
    jpype.startJVM("--enable-native-access=ALL-UNNAMED", classpath=[jar_path])
from edu.sc.seis.TauP import TauP_Time  # type: ignore


def taup_time_java(
    event_depth_km, dist_km, phases_list, receiver_depth_km=0, model_name="ak135"
):
    ttobj = TauP_Time(model_name)
    ttobj.setSourceDepth(event_depth_km)
    ttobj.setReceiverDepth(receiver_depth_km)
    ttobj.setPhaseNames(phases_list)
    ttobj.calculate(dist_km / 111.19492664455874)

    N_arr = ttobj.getNumArrivals()
    results = {"phase": [], "puristphase": [], "time": [], "rayparameter": []}
    for i in range(N_arr):
        arr = ttobj.getArrival(i)
        results["phase"].append(str(arr.getName()))
        results["puristphase"].append(str(arr.getPuristName()))
        results["time"].append(float(arr.getTime()))
        results["rayparameter"].append(float(arr.getRayParam()))
    return results


def cal_first_p(event_depth_km, dist_km, receiver_depth_km=0, model_name="ak135"):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    phases_list_p = ["p", "P", "pP", "Pg", "Pn", "Pdiff", "PKP"]
    ttobj_p = TauP_Time(model_name)
    ttobj_p.setSourceDepth(event_depth_km)
    ttobj_p.setReceiverDepth(receiver_depth_km)
    ttobj_p.setPhaseNames(phases_list_p)
    ttobj_p.calculate(dist_km / 111.19492664455874)
    arr_p = ttobj_p.getArrival(0)
    first_p = arr_p.getTime()
    return first_p


def cal_first_p_s(event_depth_km, dist_km, receiver_depth_km=0, model_name="ak135"):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    phases_list_p = ["p", "P", "pP", "Pg", "Pn", "Pdiff", "PKP"]
    ttobj_p = TauP_Time(model_name)
    ttobj_p.setSourceDepth(event_depth_km)
    ttobj_p.setReceiverDepth(receiver_depth_km)
    ttobj_p.setPhaseNames(phases_list_p)
    ttobj_p.calculate(dist_km / 111.19492664455874)
    arr_p = ttobj_p.getArrival(0)
    first_p = arr_p.getTime()

    phases_list_s = ["s", "S", "sS", "pS", "Sg", "Sn", "Sdiff", "SKS"]
    ttobj_s = TauP_Time(model_name)
    ttobj_s.setSourceDepth(event_depth_km)
    ttobj_s.setReceiverDepth(receiver_depth_km)
    ttobj_s.setPhaseNames(phases_list_s)
    ttobj_s.calculate(dist_km / 111.19492664455874)
    arr_s = ttobj_s.getArrival(0)
    first_s = arr_s.getTime()

    return first_p, first_s


def create_tpts_table(
    path_green,
    event_depth_km,
    receiver_depth_km,
    dist_km_list,
    model_name="ak135",
    check_finished=False,
):
    path_tp_table = os.path.join(
        path_green,
        "%.2f" % event_depth_km,
        "%.2f" % receiver_depth_km,
        "tp_table.bin",
    )
    path_ts_table = os.path.join(
        path_green,
        "%.2f" % event_depth_km,
        "%.2f" % receiver_depth_km,
        "ts_table.bin",
    )
    if (
        check_finished
        and os.path.exists(path_tp_table)
        and os.path.exists(path_ts_table)
    ):
        return
    tp_table = np.zeros(len(dist_km_list), dtype=np.float32)
    ts_table = np.zeros(len(dist_km_list), dtype=np.float32)
    for i in range(len(dist_km_list)):
        first_p, first_s = cal_first_p_s(
            event_depth_km, dist_km_list[i], receiver_depth_km, model_name
        )
        tp_table[i] = first_p
        ts_table[i] = first_s
    tp_table.tofile(path_tp_table)
    ts_table.tofile(path_ts_table)


if __name__ == "__main__":
    pass
