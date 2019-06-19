"""
Includes commands
"""

import argparse

__all__ = ["convert",
           "dp",
           "equal",
           "events",
           "search",
           "view_recon",
           "landscape",
           "view_landscape"]

class CustomHelpFormatter(argparse.HelpFormatter):
    """Help message formatter that ...
        adds default values to argument help.
        combines metavar for short and long options.
    """

    def _format_action_invocation(self, action):
        # See HelpFormatter._format_action_invocation

        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar

        else:
            parts = []

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s, --long ARGS
            else:
                default = action.dest.upper()
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append('%s' % option_string) # MODIFIED
                parts[-1] += ' %s' % args_string       # MODIFIED
        
            return ', '.join(parts)

    def _get_help_string(self, action):
        # ArgumentDefaultsHelpFormatter._get_help_string

        help = action.help
        if '%(default)' not in action.help:
            if (action.default is None) or \
               (action.default is True) or (action.default is False): # MODIFIED
                return help

            if action.choices is not None: # MODIFIED
                return help

            elif action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (default: %(default)s)'
                return help