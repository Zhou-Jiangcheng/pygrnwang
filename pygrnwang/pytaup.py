import os

import jpype.imports  # use jpype to call java class

jar_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "exec",
    "TauP-2.6.1.jar",
)
if not jpype.isJVMStarted():
    jpype.startJVM(classpath=[jar_path])
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
    phases_list_p = ["p", "P", "Pg", "Pn", "Pdiff", "PKP"]
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
    phases_list_p = ["p", "P", "Pg", "Pn", "Pdiff", "PKP"]
    ttobj_p = TauP_Time(model_name)
    ttobj_p.setSourceDepth(event_depth_km)
    ttobj_p.setReceiverDepth(receiver_depth_km)
    ttobj_p.setPhaseNames(phases_list_p)
    ttobj_p.calculate(dist_km / 111.19492664455874)
    arr_p = ttobj_p.getArrival(0)
    first_p = arr_p.getTime()

    phases_list_s = ["s", "S", "Sg", "Sn", "Sdiff", "SKS"]
    ttobj_s = TauP_Time(model_name)
    ttobj_s.setSourceDepth(event_depth_km)
    ttobj_s.setReceiverDepth(receiver_depth_km)
    ttobj_s.setPhaseNames(phases_list_s)
    ttobj_s.calculate(dist_km / 111.19492664455874)
    arr_s = ttobj_s.getArrival(0)
    first_s = arr_s.getTime()

    return first_p, first_s


if __name__ == "__main__":
    pass
