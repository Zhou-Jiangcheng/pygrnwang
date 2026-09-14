# Travel times and TauP

JPype is not used. Java calls execute in subprocesses; the default functions select ObsPy when the JAR or JDK tools are absent. `taup_time_java` explicitly requires Java. Importing the module detects availability without starting Java or compiling the bridge.

Java `rayparameter` is **seconds/radian**. Native SPGRN table slowness is **seconds/metre**; neither value can be used as seconds/degree without conversion.

## pytaup

```{eval-rst}
.. autofunction:: pygrnwang.pytaup.cal_first_p
```

```{eval-rst}
.. autofunction:: pygrnwang.pytaup.cal_first_s
```

```{eval-rst}
.. autofunction:: pygrnwang.pytaup.cal_first_p_s
```

```{eval-rst}
.. autofunction:: pygrnwang.pytaup.taup_time_java
```

```{eval-rst}
.. autofunction:: pygrnwang.pytaup.taup_create_npz_file
```

```{eval-rst}
.. autofunction:: pygrnwang.pytaup.create_tpts_table
```
