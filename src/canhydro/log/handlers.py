import logging
from rich.console import Console
from rich.highlighter import Highlighter, NullHighlighter
from rich.theme import Theme
from typing_extensions import Self
from typing import Any, Dict, List, Type, Union

class ConsoleHandler(logging.StreamHandler):
    def __init__(
        self,
        stream=None,
        # highlighter: Highlighter = PrefectConsoleHighlighter,
        styles: Dict[str, str] = None,
        level: Union[int, str] = logging.NOTSET,
    ):
        """
        The default console handler for Prefect, which highlights log levels,
        web and file URLs, flow and task (run) names, and state types in the
        local console (terminal).

        Highlighting can be toggled on/off with the PREFECT_LOGGING_COLORS setting.
        For finer control, use logging.yml to add or remove styles, and/or
        adjust colors.
        """
        super().__init__(stream=stream)

        if 1==1:
            highlighter = Highlighter()
            theme = Theme(styles, inherit=False)
        else:
            highlighter = NullHighlighter()
            theme = Theme(inherit=False)

        self.level = level
        self.console = Console(
            highlighter=highlighter,
            theme=theme,
            file=self.stream,
            # markup=markup_console,
        )

    def emit(self, record: logging.LogRecord):
        try:
            message = self.format(record)
            self.console.print(message, soft_wrap=True)
        except RecursionError:
            # This was copied over from logging.StreamHandler().emit()
            # https://bugs.python.org/issue36272
            raise
        except Exception:
            self.handleError(record)