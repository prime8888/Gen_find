# function.py
from PyQt5.QtWidgets import QCheckBox
from PyQt5.QtCore import Qt

def toggle_all_checkboxes(checkboxes, state, log_text):
    """Toggle all given checkboxes according to the state."""
    for checkbox in checkboxes:
        checkbox.blockSignals(True)
        checkbox.setChecked(state == Qt.Checked)
        checkbox.blockSignals(False)
    if state == Qt.Checked:
        log_text.append("All checkboxes selected.")
    else:
        log_text.append("All checkboxes deselected.")

def checkbox_toggled(option, state, checkbox, all_checkbox, log_text):
    """Handle the toggling of an individual checkbox and adjust 'ALL' checkbox accordingly."""
    toggle_status = 'checked' if state else 'unchecked'
    log_text.append(f"{option} checkbox toggled to {toggle_status}.")

    if not state:
        all_checkbox.blockSignals(True)
        all_checkbox.setChecked(False)
        all_checkbox.blockSignals(False)
    else:
        # Check if all checkboxes (except 'ALL') are checked
        checkboxes = checkbox.parent().findChildren(QCheckBox)
        if all(chk.isChecked() for chk in checkboxes if chk != all_checkbox):
            all_checkbox.blockSignals(True)
            all_checkbox.setChecked(True)
            all_checkbox.blockSignals(False)
