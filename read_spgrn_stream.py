from obspy import Trace, Stream, UTCDateTime
from pygrnwang.read_spgrn_old import Spgrn
from pygrnwang.signal_process import filter_butter, resample


def read_spgrn_obspy_stream(
    path_green_lib,
    event_time="",
    event_depth_in_km=None,
    az_in_deg=None,
    dist_in_km=None,
    focal_mechanism=None,
    before_tP=120,
    time_window=1024,
    srate=None,
    freq_min=0,
    freq_max=0,
    tension_compression_flag=False,
):
    spgrn = Spgrn(
        path_green_lib=path_green_lib,
        event_depth_in_km=event_depth_in_km,
        az_in_deg=az_in_deg,
        dist_in_km=dist_in_km,
        focal_mechanism=focal_mechanism,
        before_tp=before_tP,
        time_window=time_window,
        tension_compression_flag=tension_compression_flag,
    )
    spgrn.read_enz()
    stream = Stream()
    for i in range(3):
        trace = Trace(data=spgrn.waveforms[i])
        trace.stats.sampling_rate = 1 / spgrn.green_info["sampling_interval"]
        if event_time != "":
            event_time = UTCDateTime(event_time)
            trace.stats.starttime = event_time + spgrn.green_tpts_table["p_onset"]

        trace.data = filter_butter(
            data=trace.data,
            srate=1 / spgrn.green_info["sampling_interval"],
            freq_band=[freq_min, freq_max],
            butter_order=4,
            zerophase=False,
        )
        trace.resample(sampling_rate=srate)
        stream.append(trace)
    return stream


if __name__ == "__main__":
    path = "/e/spgrn2020_lib"
    az_ = 0
    dist_in_km_ = 5000
    fm_ = [30, 50, 60]
    st = read_spgrn_obspy_stream(
        path_green_lib=path,
        event_time="2020-03-24T20:00:12.123",
        event_depth_in_km=10,
        az_in_deg=az_,
        dist_in_km=dist_in_km_,
        focal_mechanism=fm_,
        before_tP=120,
        time_window=1024,
        srate=4,
        freq_min=0,
        freq_max=1,
        tension_compression_flag=False,
    )
    st[-1].plot()
