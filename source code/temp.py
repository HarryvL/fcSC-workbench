from PySide2.QtWidgets import QApplication, QMainWindow, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget, \
    QLineEdit, QPushButton
from PySide2.QtCore import Qt
import sys


class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()

        self.setWindowTitle("Dynamic QTableWidget Example with Editable and Non-Editable Cells")
        self.setGeometry(100, 100, 600, 400)

        # Create a QTableWidget
        self.table_widget = QTableWidget()

        # Set initial column count
        self.table_widget.setColumnCount(3)

        # Set column headers
        self.table_widget.setHorizontalHeaderLabels(['Column 1', 'Column 2', 'Column 3'])

        # Layout to hold the table and button
        layout = QVBoxLayout()
        layout.addWidget(self.table_widget)

        # Create a button to add new rows
        self.add_row_button = QPushButton("Add Row")
        self.add_row_button.clicked.connect(self.add_row)
        layout.addWidget(self.add_row_button)

        # Create a central widget and set the layout
        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

        # Dynamically populate the table
        self.populate_table()

    def populate_table(self):
        data = [
            ["Label1", "Input1", "Label2"],
            ["Label3", "Input2", "Label4"],
            ["Label5", "Input3", "Label6"]
        ]

        for row_data in data:
            self.add_row(row_data)

    def add_row(self, row_data=None):
        if row_data is None:
            # Default data for new row
            row_data = ["New Label", "New Input", "New Label"]

        row_index = self.table_widget.rowCount()
        self.table_widget.insertRow(row_index)

        for column_index, item in enumerate(row_data):
            if column_index % 2 == 0:
                # Non-editable cell (label)
                table_item = QTableWidgetItem(item)
                table_item.setFlags(table_item.flags() & ~Qt.ItemIsEditable)
                self.table_widget.setItem(row_index, column_index, table_item)
            else:
                # Editable cell (user input)
                line_edit = QLineEdit()
                line_edit.setText(item)
                self.table_widget.setCellWidget(row_index, column_index, line_edit)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
