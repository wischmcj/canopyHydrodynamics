import logging 
import yaml
from rich.console import Console
from rich.highlighter import Highlighter
from rich.theme import Theme
from typing import Union


with open('./src/canhydro/logging_config.yml', 'rt') as f:
    config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

class ConsoleHandler(logging.StreamHandler):
    def __init__(
        self,
        stream=None,
        highlighter: Highlighter = Highlighter,
        styles: dict[str, str] = None,
        level: Union[int, str] = logging.NOTSET,
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

        highlighter = highlighter()
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
            print(f'Error in console handler: {e}')