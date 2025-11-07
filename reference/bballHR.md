# MLB Players' Home Run and Batted Ball Statistics with Red Zone Metrics (2017-2019)

This dataset provides annual statistics for Major League Baseball (MLB)
players, including home run counts, at-bats, mean exit velocities,
launch angles, quantile statistics of exit velocities and launch angles,
and red zone metrics. It is intended for analyzing batted ball
performance, with additional variables on the red zone, which is defined
as balls in play with a launch angle between 20 and 35 degrees and an
exit velocity of at least 95 mph.

## Usage

``` r
bballHR
```

## Format

A data frame with the following columns:

- name:

  Player's full name (character).

- playerID:

  Player's unique identifier in the Lahman database (character).

- teamID:

  Team abbreviation (character).

- year:

  Season year (numeric).

- HR:

  Home runs hit during the season (integer).

- AB:

  At-bats during the season (integer).

- mean_exit_velo:

  Average exit velocity (mph) over the season (numeric).

- mean_launch_angle:

  Average launch angle (degrees) over the season (numeric).

- launch_angle_75:

  Launch angle at the 75th percentile of the player's distribution
  (numeric).

- launch_angle_70:

  Launch angle at the 70th percentile of the player's distribution
  (numeric).

- launch_angle_65:

  Launch angle at the 65th percentile of the player's distribution
  (numeric).

- exit_velo_75:

  Exit velocity at the 75th percentile of the player's distribution
  (numeric).

- exit_velo_80:

  Exit velocity at the 80th percentile of the player's distribution
  (numeric).

- exit_velo_85:

  Exit velocity at the 85th percentile of the player's distribution
  (numeric).

- count_red_zone:

  Seasonal count of batted balls in the red zone, defined as a launch
  angle between 20 and 35 degrees and an exit velocity greater than or
  equal to 95 mph (integer).

- prop_red_zone:

  Proportion of batted balls that fall into the red zone (numeric).

- BPF:

  Ballpark factor, indicating the effect of the player's home ballpark
  on offensive statistics (integer).

## Source

- Player statistics: [Lahman R
  Package](https://CRAN.R-project.org/package=Lahman)

- Batted ball data: [Baseball Savant](https://baseballsavant.mlb.com/)

- Additional analysis: [Patterns of Home Run Hitting in the Statcast
  Era](https://bayesball.github.io/BLOG/homeruns.html) by Jim Albert

## Details

- Mean Metrics:

  `mean_exit_velo` and `mean_launch_angle` represent the player's
  average exit velocities and launch angles, respectively, over the
  course of a season.

- Quantile Metrics:

  The `launch_angle_xx` and `exit_velo_xx` columns denote the upper
  \\x\\-percentiles (e.g., 75th percentile) of the player's launch angle
  and exit velocity distributions for that year.

- Red Zone Metrics:

  `count_red_zone` gives the number of balls in play that fall into the
  red zone, while `prop_red_zone` represents the proportion of balls in
  play in this category.

- BPF:

  The Ballpark Factor (BPF) quantifies the influence of the player's
  home ballpark on offensive performance, with values above 100
  indicating a hitter-friendly environment.

## Examples

``` r
data(bballHR)
head(bballHR)
#>            name  playerID teamID year HR  AB mean_exit_velo mean_launch_angle
#> 1    Jose Abreu abreujo02    CHA 2017 33 621       90.54188          11.27789
#> 2    Jose Abreu abreujo02    CHA 2018 22 499       91.31206          12.08543
#> 3    Jose Abreu abreujo02    CHA 2019 33 634       92.12946          10.97996
#> 4  Ronald Acuna acunaro01    ATL 2019 41 626       90.65088          14.04825
#> 5  Willy Adames adamewi01    TBA 2019 20 531       88.38465          10.50128
#> 6 Jesus Aguilar aguilje01    MIL 2018 35 492       89.78601          16.29793
#>   launch_angle_75 launch_angle_70 launch_angle_65 exit_velo_75 exit_velo_80
#> 1              27            23.0           20.00      102.600       104.30
#> 2              27            24.0           20.00      102.900       104.32
#> 3              26            22.0           18.00      103.000       104.30
#> 4              30            28.0           26.00      103.525       104.90
#> 5              27            23.0           20.00       99.350       101.00
#> 6              33            30.5           26.25      100.375       102.10
#>   exit_velo_85 count_red_zone prop_red_zone BPF
#> 1      105.750             74     0.1312057  98
#> 2      106.135             53     0.1247059  97
#> 3      106.000             76     0.1404806  97
#> 4      106.675             86     0.1787942 105
#> 5      102.650             62     0.1486811  97
#> 6      103.325             71     0.1674528 102
```
