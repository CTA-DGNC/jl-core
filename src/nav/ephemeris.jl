module ephemeris

#=
  For the Model E AR5 simulations, all years have exactly 365 days and 
  vernal equinox always occurs on March 21, hour 0. 
  Future versions of Model E will use the more precise time of vernal 
  equinox as described on this web page. 

  The date of vernal equinox and the number of days in February are 
  political decisions when choosing to use the Gregorian calendar. 
  In this calendar, February has 28 days most years every fourth year 
  (1988 A.D., 1992, 1996) February has 29 days every century (1700, 
  1800, 1900) February has 28 days and every fourth century (1600, 
  2000, 2400) February has 29 days. ephemeris

  This web page assumes that this 400 year cycle is repeated indefinitely, 
  that the tropical year is (365*400+97)/400 = 365.2425 days, and that 
  vernal equinox occurs exactly on March 20, 7:30 GMT every four 
  hundred years (including year 2000 A.D.). 

  This last decision is based on comparing the web page's time of vernal 
  equinox with observations. The Table shows web page times and actual 
  times of vernal equinox published in NASA Reference Publication 1349 
  [1994 October] for years 1995 to 2010 and in Explanatory Statement 
  to the Ephemeris [19??] for years 1903, 2000 and 2096.
=#

    const Jan_01_1601 = 584754
    const Jan_01_1900 = 693961
    const Jan_01_1970 = 719528

    const month_days       = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    const days_up_to_month = [ 0
                            31
                            31 + 28
                            31 + 28 + 31
                            31 + 28 + 31 + 30
                            31 + 28 + 31 + 30 + 31
                            31 + 28 + 31 + 30 + 31 + 30
                            31 + 28 + 31 + 30 + 31 + 30 + 31
                            31 + 28 + 31 + 30 + 31 + 30 + 31 + 31
                            31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30
                            31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31
                            31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30]

    # https://data.giss.nasa.gov/modelE/ar5plots/srvernal.html
    const ♈_utc = [2000 3 20 7 30 0] 
    const ♈_day = 730199
    const ♈_sec = ♈_day * 86400

    
    function leap_count(nYear) 
        nYear = nYear - 1
        floor(nYear / 4) - floor(nYear / 100) + floor(nYear / 400)
    end
    function is_leap_year(nYear)  
        ~mod(nYear, 4) && (mod(nYear, 100) || ~mod(nYear, 400))
    end
    function get_day_count(nYear, eMonth, nDay)
        nDayCount = (nYear - 1) * 365 + leap_count(nYear) + nDay
        nDayCount = nDayCount + days_up_to_month[eMonth]
        if eMonth > 2 && is_leap_year(nYear)  
            nDayCount = nDayCount + 1
        end
        nDayCount
    end
    # segundos o [y m d hh mm ss]
    function UT(utc)
        if size(utc) > 2
            tm = get_day_count(utc[1], utc[2], utc[3]) * 86400
            if length(utc) > 3
                tm = tm + utc[4] * 3600 
            end
            if length(utc) > 4
                tm = tm + utc[5] * 60 
            end
            if length(utc) > 6
                tm = tm + utc[6] 
            end
        else
            tm = utc
        end
        tm - ♈_sec
    end        

end
