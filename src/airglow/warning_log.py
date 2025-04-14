from datetime import datetime

class WarningLog:
    def __init__(self):
        self.log = []
        self.warning_labels = set()    # Using set automatically handles duplicates
        self.multiple_warnings = False  # Flag to track multiple different warnings

    def add(self, message: str, warning_type: str = None, title: str = None, label: str = None):
        timestamp = datetime.now().strftime('%m/%d/%Y %H:%M:%S %p')
        self.log.append(f"{timestamp}: {message}")
        
        # Track if this is a subsequent warning with a different title
        if title and hasattr(self, 'first_title'):
            if title != self.first_title:
                self.multiple_warnings = True
        else:
            self.first_title = title
        
        # Add label if provided and not already present
        if label:
            self.warning_labels.add(label)  # set.add() automatically handles duplicates

    @property
    def warning_title(self):
        """Return appropriate title based on number of different warnings"""
        if self.multiple_warnings:
            return "Multiple warnings encountered"
        return getattr(self, 'first_title', None)

    def get_full_log(self) -> str:
        return "\n".join(self.log)
    
    def __str__(self) -> str:
        return self.get_full_log()
    
    def __bool__(self) -> bool:
        return len(self.log) > 0