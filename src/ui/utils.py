# function.py
from PyQt5.QtWidgets import QCheckBox
from PyQt5.QtCore import Qt

def toggle_all_checkboxes(checkboxes, state):
    """Toggle all given checkboxes according to the state."""
    for checkbox in checkboxes:
        checkbox.blockSignals(True)
        checkbox.setChecked(state == Qt.Checked)
        checkbox.blockSignals(False)

def checkbox_toggled(option, state, checkbox, all_checkbox):
    """Handle the toggling of an individual checkbox and adjust 'ALL' checkbox accordingly."""
    toggle_status = 'checked' if state else 'unchecked'

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
