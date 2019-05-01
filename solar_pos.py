# This should calculate the solar position baswed on geographic position and time of year

import numpy as np
import pandas as pd
import os
import sys
import math
import datetime

def calc_sun_position(latitude_deg, longitude_deg, year, hoy):
        """
        Calculates the Sun Position for a specific hour and location
        :param latitude_deg: Geographical Latitude in Degrees
        :type latitude_deg: float
        :param longitude_deg: Geographical Longitude in Degrees
        :type longitude_deg: float
        :param year: year
        :type year: int
        :param hoy: Hour of the year from the start. The first hour of January is 1
        :type hoy: int
        :return: altitude, azimuth: Sun position in altitude and azimuth degrees [degrees]
        :rtype: tuple
        """

        # Convert to Radians
        latitude_rad = math.radians(latitude_deg)
        longitude_rad = math.radians(longitude_deg)

        # Set the date in UTC based off the hour of year and the year itself
        start_of_year = datetime.datetime(year, 1, 1, 0, 0, 0, 0)
        utc_datetime = start_of_year + datetime.timedelta(hours=hoy)

        # Angular distance of the sun north or south of the earths equator
        # Determine the day of the year.
        day_of_year = utc_datetime.timetuple().tm_yday

        # Calculate the declination angle: The variation due to the earths tilt
        # http://www.pveducation.org/pvcdrom/properties-of-sunlight/declination-angle
        declination_rad = math.radians(
            23.45 * math.sin((2 * math.pi / 365.0) * (day_of_year - 81)))

        # Normalise the day to 2*pi
        # There is some reason as to why it is 364 and not 365.26
        angle_of_day = (day_of_year - 81) * (2 * math.pi / 364)

        # The deviation between local standard time and true solar time
        equation_of_time = (9.87 * math.sin(2 * angle_of_day)) - \
            (7.53 * math.cos(angle_of_day)) - (1.5 * math.sin(angle_of_day))

        # True Solar Time
        solar_time = ((utc_datetime.hour * 60) + utc_datetime.minute +
                      (4 * longitude_deg) + equation_of_time) / 60.0

        # Angle between the local longitude and longitude where the sun is at
        # higher altitude
        hour_angle_rad = math.radians(15 * (12 - solar_time))

        # Altitude Position of the Sun in Radians
        altitude_rad = math.asin(math.cos(latitude_rad) * math.cos(declination_rad) * math.cos(hour_angle_rad) +
                                 math.sin(latitude_rad) * math.sin(declination_rad))

        # Azimuth Position fo the sun in radians
        azimuth_rad = math.asin(
            math.cos(declination_rad) * math.sin(hour_angle_rad) / math.cos(altitude_rad))

        # I don't really know what this code does, it has been imported from
        # PySolar
        if(math.cos(hour_angle_rad) >= (math.tan(declination_rad) / math.tan(latitude_rad))):
            return math.degrees(altitude_rad), math.degrees(azimuth_rad)
        else:
            return math.degrees(altitude_rad), (180 - math.degrees(azimuth_rad))
