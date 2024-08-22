from __future__ import annotations

import logging
from typing import ClassVar

from rich.console import Console
from rich.highlighter import Highlighter, ReprHighlighter
from rich.theme import Theme


class ConsoleHandler(logging.StreamHandler):
    HIGHLIGHTER_CLASS: ClassVar[type[Highlighter]] = ReprHighlighter

    def __init__(
        self,
        stream=None,
        highlighter: Highlighter | None = None,
        styles: dict[str, str] = None,
        level: int | str = logging.NOTSET,
    ):
        """
        Allows for logging to the console with rich formatting.
        Args:
            stream (Optional[IO[str]]): A file-like object to write to. Defaults to sys.stderr.
            highlighter (Highlighter, optional): A highlighter to use for syntax
                highlighting. Using  NullHighlighter eliminates highlighting.
            styles (Dict[str, str], optional): A dictionary of styles to use. Defaults to None.
            level (Union[int, str], optional): The minimum log level to write to the console.
                Defaults to logging.NOTSET.

        Note: 'markup = True' allows for customization of error message format based on bracketed tags
                This may cause strange behavior with brackets
        """
        super().__init__(stream=stream)

        highlighter = highlighter or self.HIGHLIGHTER_CLASS()
        theme = Theme(styles, inherit=False)

        self.level = level
        self.console = Console(
            highlighter=highlighter,
            theme=theme,
            file=self.stream,
            markup=True,
        )

    def emit(self, record: logging.LogRecord):
        try:
            message = self.format(record)
            self.console.print(message, soft_wrap=True)
        except Exception as e:
            print(f"Error in console handler: {e}")
