"""
A script to sign in to Studiosity and check links
"""
import sys
import time
import random
from selenium import webdriver
from selenium.common.exceptions import NoSuchElementException

email_login = ""
password = ""

# Open the driver with Chrome
chrome_driver = "/Users/Chow/Documents/computational/Python/chromedriver"
driver = webdriver.Chrome(chrome_driver)

# Try to go straight to tutor page incase already signed in, otherwise login
driver.get("https://studiosity.com/connect/tutors")
if driver.title == "Sign In - Studiosity":
    if email_login and password:
        # Select the email login box and send key
        email_box = driver.find_element_by_id("user_email")
        email_box.send_keys(email_login)
        time.sleep(random.randint(1,3))

        # Select password and send key
        pass_box = driver.find_element_by_id("user_password")
        pass_box.send_keys(password)
        time.sleep(random.randint(1,3))

        # Login button
        login_button = driver.find_element_by_name("commit")
        login_button.click()
    else:
        hold = input("Enter your login details, sign in, then press Enter")

# should be on page studiosity.com/connect/tutors
assert 'Specialist Home' in driver.title

# check the submission button, if a submission is current then this will exit
i = 0
while True:
    try:
        sub_box = driver.find_element_by_partial_link_text(
                                            "View Writing Feedback Home - ")
    except NoSuchElementException:
        try:
            driver.find_element_by_partial_link_text("Current Submission - ")
            print("Current Submission, please try again later.")
            sys.exit()
        except NoSuchElementException:
            print("Error 1: Couldn't find submission button")

    # check the submission count 
    try:
        count_elem = sub_box.find_element_by_class_name(
                                            "submission-count-description")
    except NoSuchElementException:
        print("Error 2: Couldn't find submission count description")

    # get the description
    submission_count = int(count_elem.text.split()[0])

    # if there are 0 submissions, sleep some random amount of time, else click
    if submission_count == 0:
        random_number = random.randint(1,6)
        print(f"0 Submissions, waiting {random_number}s to refresh")
        time.sleep(random_number)
        driver.refresh()
        i += 1
        if i == 25:
            check = print("Do you want to continue refreshing? ([y]/n): ")
            if check == "n":
                sys.exit()
            else:
                i = 0
    else:
        print("click")
        sub_box.click()
        hold = input("Submission found, do not close this window!")
        hold = input("Don't push enter!")
        hold = input("I mean it!")
        sys.exit()






