"""Headless browser smoke test for index.html (offline, file://).
Run: python test/smoke_browser.py
Checks: page loads with no severe console errors, each tier reconstructs, badge renders,
export is enabled for A/B and refused for C.
"""
import io
import os
import sys
import time
import pathlib

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC

HTML = pathlib.Path(__file__).resolve().parent.parent / "index.html"
URL = HTML.as_uri()


def main():
    opts = Options()
    opts.add_argument("--headless=new")
    opts.add_argument("--allow-file-access-from-files")
    opts.add_argument("--disable-gpu")
    opts.set_capability("goog:loggingPrefs", {"browser": "ALL"})
    driver = webdriver.Chrome(options=opts)
    failures = []
    try:
        driver.set_page_load_timeout(60)
        driver.get(URL)
        wait = WebDriverWait(driver, 30)
        # auto-loads first example into the paste box
        wait.until(lambda d: d.find_element(By.ID, "paste").get_attribute("value").strip() != "")

        def reconstruct(example_key):
            Select(driver.find_element(By.ID, "example")).select_by_value(example_key)
            driver.find_element(By.ID, "loadEx").click()
            driver.find_element(By.ID, "run").click()
            wait.until(EC.visibility_of_element_located((By.ID, "out")))
            time.sleep(0.4)
            tier = driver.find_element(By.ID, "tier").text
            badge = driver.find_element(By.ID, "badge").text
            export_disabled = driver.find_element(By.ID, "exportBtn").get_attribute("disabled")
            checks = driver.find_elements(By.CSS_SELECTOR, "#checks tbody tr")
            paths = driver.find_elements(By.CSS_SELECTOR, "#plot path")
            return tier, badge, export_disabled, len(checks), len(paths)

    # Tier A
        tier, badge, exp_dis, nchecks, npaths = reconstruct("tierA")
        print(f"Tier A -> {tier}, badge={badge}, exportDisabled={exp_dis}, checks={nchecks}, curves={npaths}")
        if "A" not in tier: failures.append("Tier A not classified A")
        if exp_dis: failures.append("Tier A export should be enabled")
        if npaths < 2: failures.append("Tier A should draw 2 reconstructed curves")
        if nchecks < 5: failures.append("Tier A checks table missing rows")

        # Tier B
        tier, badge, exp_dis, nchecks, npaths = reconstruct("tierB")
        print(f"Tier B -> {tier}, badge={badge}, exportDisabled={exp_dis}, checks={nchecks}, curves={npaths}")
        if "B" not in tier: failures.append("Tier B not classified B")
        if exp_dis: failures.append("Tier B export should be enabled")

        # Tier C
        tier, badge, exp_dis, nchecks, npaths = reconstruct("tierC")
        print(f"Tier C -> {tier}, badge={badge}, exportDisabled={exp_dis}")
        if "C" not in tier: failures.append("Tier C not classified C")
        if not exp_dis: failures.append("Tier C export MUST be refused (disabled)")
        if badge.lower() != "none": failures.append("Tier C badge must be NONE")

        # console errors (ignore favicon)
        severe = [e for e in driver.get_log("browser")
                  if e["level"] == "SEVERE" and "favicon" not in e["message"]]
        for e in severe:
            failures.append("console SEVERE: " + e["message"][:200])
    finally:
        driver.quit()

    if failures:
        print("\nFAIL:")
        for f in failures:
            print("  -", f)
        return 1
    print("\nALL BROWSER SMOKE CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
