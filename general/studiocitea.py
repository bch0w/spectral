"""
A script to sign in to Studiosity and check links
"""
import os
import sys
import time
import random
import getpass
from selenium import webdriver
from selenium.common.exceptions import NoSuchElementException


def login(driver, email_login='', password=''):
    """
    Login to the webpage either with a given email/password, or prompted to
    enter credentials yourself

    :type driver: selenium.webdriver.chrome.webdriver.WebDriver
    :param driver: the webdriver that is accessing the page
    :type email_login: str
    :type email_login: email login for webpage, requires @. If blank, function
        will prompt a manual login
    :type password: str
    :param password: plain text password
    """
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
    if 'Specialist Home' in driver.title:
        return 1
    else:
        print("incorrect password")
        return 0


def refresh_home(driver, refresh_max=25, sleep_range=(1,7)):
    """
    Sit on the home page and refresh the page and check a certain element to
    determine the number of submissions, if it goes above a certain number
    then click the submission button

    TO DO:
        Can the element finding be moved outside the loop?

    :type driver: selenium.webdriver.chrome.webdriver.WebDriver
    :param driver: chrome driver for controlling web browser
    :type refresh_max: int
    :param refresh_max: max refresh limit before prompting response
    :type sleep_range: tuple of int
    :param sleep_range: range for the R.N.G to choose a sleep time
    """
    # make sure were on the right page
    if 'Specialist Home' not in driver.title:
        driver.get("https://studiosity.com/connect/tutors")

    refreshed = 0
    while True:
        # find the button to check
        try:
            sub_box = driver.find_element_by_partial_link_text(
                                                "View Writing Feedback Home - ")
        # button may not be available if a submission is current
        except NoSuchElementException:
            try:
                driver.find_element_by_partial_link_text("Current Submission -")
                print("Current Submission, please try again later.")
                return -1
            except NoSuchElementException:
                print("Error 1: Couldn't find submission button")
                return -1

        # check the submission count 
        try:
            count_elem = sub_box.find_element_by_class_name(
                                                "submission-count-description")
        except NoSuchElementException:
            print("Error 2: Couldn't find submission count description")
            return -1

        # get the submission count as an int
        submission_count = int(count_elem.text.split()[0])

        # if 0 submissions, sleep some random amount of time, or click
        if submission_count == 0:
            # check how many times the page has been refreshed
            refreshed += 1
            if refreshed == refresh_max:
                return 0

            random_number = random.randint(sleep_range[0], sleep_range[1])
            print(f"0 submissions - waiting {random_number}s to refresh")
            time.sleep(random_number)
            driver.refresh()
        else:
            print("click")
            # sub_box.click()
            os.system("say -v Thomas 'submission found'")
            # TO DO, continue to find elements here
            return 1
    

def check():
    """
    Since closing the script shuts the webpage, this function will hold the 
    script open while work is being done
    """
    check_ = input("Waiting for response...: ")
    return check_


def state_machine(driver, refresh_max=25):
    """
    The main workflow that allows indefinite use of this script during work

    :type driver: selenium.webdriver.chrome.webdriver.WebDriver
    :param driver: chrome driver for controlling web browser
    :type refresh_max: int
    :param refresh_max: max refresh limit before prompting response
    """
    def submit_state():
        """
        State defined by when a submission has been found
        """
        state_submit = None
        while state_submit not in ["y", "n"]:
            state_submit = input("Submission found, accepted? y/n: ")
        if state == "y":
            print("Once submission completed, press any key to resume")
            check()
        elif state == "n":
            print("Continuing to refresh")
        return 1

    def continue_state():
        """
        State where User intervention is required to continue
        """
        state_continue = None 
        while state_continue not in ["1", "2"]:
            state_continue = input("What would you like to do?:\n"
                                   "\t1) Exit\n"
                                   "\t2) Continue refreshing\n")
        if state_continue == "1":
            sys.exit("Goodbye")
        else:
            print("Continuing to refresh")
        return 1

    def refresh_state():
        """
        State the determine if refresh should be continued
        """
        state_refresh = None
        while state_refresh not in ["y", "n"]:
            state_refresh = input("Continue refreshing? y/n: ")
        if state == "y":
            print("Continuing to refresh")
            return 1
        elif state == "n":
            return 0

    try:
        # Open the driver with Chrome, login
        email_login = "edwardsamj@live.com"

        login_success = False
        while not login_success:
            password = getpass.getpass(f"Password for '{email_login}' please: ")
            login_status = login(driver, email_login, password)
            if login_status:
                login_success = True
        print("\n")
        
        while True:
            submission = refresh_home(driver, refresh_max)
            # Submission found
            if submission == 1:
                submit_state()
            # Refresh max reached
            elif submission == 0:
                refresh = refresh_state()
                if not refresh:
                    continue_state()
            # Some error condition has been hit
            elif submission == -1:
                continue_state()
    except KeyboardInterrupt:
        sys.exit("\nGoodbye")
    

if __name__ == "__main__":
    refresh_max = 25
    chrome_driver = webdriver.Chrome("/Users/Chow/Documents/computational/Python/chromedriver")
    state_machine(chrome_driver, refresh_max)

