 # Elephants-Moving-Data
A notebook for exploratory data analysis

First, a few graphs are linked below to get an overview of the data set: <br> <br>
[Overview with an ordinary map](https://timom2110.github.io/Elephants-Moving-Data/map_chart.html)

[Overview with a satellite view](https://timom2110.github.io/Elephants-Moving-Data/map_chart2.html)

[Yet another overview with Satellite and Streets View](https://timom2110.github.io/Elephants-Moving-Data/map_chart3.html)

---
in case Github won't show the Jupyter Notebook (.ipynb):
https://nbviewer.org/github/TimoM2110/Elephants-Moving-Data/blob/master/elephantsexploration.ipynb

---

A few plots that emerged from the exploratory data analysis (see .ipynb scripts)

---

The most important things first: This is an unevenly spaced Dataset, as you can see in [this](https://timom2110.github.io/Elephants-Moving-Data/intervals_between_records.html) as well as in [this](https://timom2110.github.io/Elephants-Moving-Data/number_of_records_per_day.html) plot. <br>

[Also the amount of records differ between the 14 tracked herds.](https://timom2110.github.io/Elephants-Moving-Data/countplot.svg) <br>

[Speed distribution](https://timom2110.github.io/Elephants-Moving-Data/Histogram_of_speeds.html) and the [directions](https://timom2110.github.io/Elephants-Moving-Data/Histogram_of_directions.html) do not show any irregular pattern though. <br>

What is more interesting is the [mean speed by temperature and month of year](https://timom2110.github.io/Elephants-Moving-Data/Mean_speed_by_temperature_and_month_of_year.html)

A HMM will be fitted in order to examine different states.
