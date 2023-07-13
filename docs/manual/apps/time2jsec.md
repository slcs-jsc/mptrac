# time2jsec

Given a date and time this app returns the Julian seconds since 01.01.2000.

```
# calling time2jsec
$ ./time2jsec  <year> <month> <day> <hour> <minute> <second> <remain>
```
The required arguments are:
* year: The year, e.g. 2020.
* month: Valid parameters for month are 1–12.
* day: Valid values for day are 1–31.
* hour: Valid values for hour are 0–23.*
* minute: Valid values for minute are 0–59.
* second: Valid values for second are 0–59.
* remain: Can be used add milliseconds, e.g. 0.111. If this value is set to 1, a second is added.

**Note:** This app also accepts values outside the valid value range whithout any warning. 
